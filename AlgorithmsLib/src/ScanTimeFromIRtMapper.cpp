//
// Created by anichols on 8/24/23.
//

#include "ScanTimeFromIRtMapper.h"

#include "CalibrationReader.h"
#include "ErrorUtils.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

ScanTimeFromIRtMapper::ScanTimeFromIRtMapper()
: m_segments(0)
, m_rtSegments(20)
, m_minRtPredBin(20)
{}


Err ScanTimeFromIRtMapper::init(const QString &iRTRecalibrationFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(iRTRecalibrationFilePath); ree;

    QVector<IRTReCalibrationRow> iRTReCalibrationReaderRows;
    e  = CSVReader::read(
            iRTRecalibrationFilePath,
            &iRTReCalibrationReaderRows
    ); ree;

    const auto insertLogic
        = [](const IRTReCalibrationRow &r){return QPair<double, double>(r.iRT, r.scanTime);};

    QVector<QPair<double, double>> data;
    std::transform(
            iRTReCalibrationReaderRows.begin(),
            iRTReCalibrationReaderRows.end(),
            std::back_inserter(data),
            insertLogic
            );

    m_segments = std::min(
            m_rtSegments,
            static_cast<int>(std::max(1.0, 2.0 * sqrt(data.size() / m_minRtPredBin)))
            );

    e = mapRT(data); ree;

    ERR_RETURN
}


namespace {

    // implementation of the algorithm from
    // stat.wikia.com/wiki/Isotonic_regression
    Err PAVA(
            const QVector<QPair<double, double>> &data,
            QVector<QPair<double, double>> *resultOutput
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(data); ree;

        QVector<int> pavaIndex;
        QVector<double> pavaW, pavaA;

        int i;
        int j;
        int k;
        int l;
        int dataSize = data.size();

        pavaIndex.resize(dataSize + 1);
        pavaW.resize(dataSize);
        pavaA.resize(dataSize);

        QVector<QPair<double, double>> result(dataSize);
        for (i = 0; i < dataSize; i++) {
            result[i].first = data[i].first;
        }

        pavaA[0] = data[0].second;
        pavaW[0] = 1.0;
        j = 0;
        pavaIndex[0] = 0;
        pavaIndex[1] = 1;

        for (i = 1; i < dataSize; i++) {

            j++;

            pavaA[j] = data[i].second;
            pavaW[j] = 1.0;
            while (j > 0 && pavaA[j] < pavaA[j - 1]) {

                const double denom = (pavaW[j] + pavaW[j - 1]);
                e = ErrorUtils::isFalse(MathUtils::tZero(denom)); ree;

                pavaA[j - 1] = (pavaW[j] * pavaA[j] + pavaW[j - 1] * pavaA[j - 1]) / denom;
                pavaW[j - 1] = pavaW[j] + pavaW[j - 1];
                j--;
            }

            pavaIndex[j + 1] = i + 1;
        }

        for (k = 0; k <= j; k++) {
            for (l = pavaIndex[k]; l < pavaIndex[k + 1]; l++) {
                result[l].second = pavaA[k];
            }
        }

        *resultOutput = result;

        ERR_RETURN
    }

    Err spline(
            const QVector<QPair<double, double>> &data,
            int segments,
            QVector<double> *coeffs,
            QVector<double> *points
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(data); ree;

        double w;

        int i;
        int k;
        int m;
        const int dataSize = data.size();

        using Triplet = Eigen::Triplet<double>;
        std::vector<Triplet> TrL;

        QVector<double> rtp;
        QVector<double> rte;
        QVector<double> rts;

        QVector<QPair<double, double>> pava;
        e = PAVA(data, &pava); ree;

        rtp.resize(dataSize);
        rte.resize(dataSize);

        m = segments;

        coeffs->resize(2 * m);

        rts.resize(m + 1);
        points->resize(m);

        for (i = 0; i < dataSize; i++) {
            rtp[i] = data[i].second, rte[i] = data[i].first;
        }

        e = ErrorUtils::isNotEqual(m, 0); ree
        const double split = (static_cast<double >(dataSize - 1)) / static_cast<double>(m);

        for (i = 0; i < m; i++) {
            rts[i] = rte[std::min(static_cast<int>(split * i), dataSize - 1)];
            rts[m] = rte[dataSize - 1];
        }

        for (i = 0; i < m; i++) {
            (*points)[i] = 0.5 * (rts[i] + rts[i + 1]);
        }

        for (i = 1; i < m; i++) {
            if (points->at(i) - points->at(i - 1) < std::numeric_limits<double>::min()) {
                for (k = i - 1; k < m - 1; k++) {
                    points[k] = points[k + 1];
                }
                m--;
            }
        }

        TrL.clear();

        for (i = k = 0; i < dataSize; i++) {
            double x = rte[i];

            if (!k && x <= points->at(i)) {
                TrL.emplace_back(i, 0, 1.0);
                TrL.emplace_back(i, 1, x - points->at(0));
            }
            else {
                for (; k < m && points->at(k) < x; k++);
                if (k == m) {
                    TrL.emplace_back(i, (k - 1) * 2, 1.0);
                    TrL.emplace_back(i, (k - 1) * 2 + 1, x - points->at(k - 1));
                }
                else {
                    w = points->at(k) - points->at(k - 1);
                    e = ErrorUtils::isFalse(MathUtils::tZero(w)); ree

                    const double u = (x - points->at(k - 1)) / w;
                    const double v = u - 1.0;

                    TrL.emplace_back(i, (k - 1) * 2, (1.0 + 2.0 * u) * std::pow(v, 2));
                    TrL.emplace_back(i, (k - 1) * 2 + 1, w * u * std::pow(v, 2));

                    TrL.emplace_back(i, k * 2, std::pow(u, 2) * (1.0 - 2.0 * v));
                    TrL.emplace_back(i, k * 2 + 1, w * std::pow(u, 2) * v);
                }
            }
        }

        const Eigen::VectorXd B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rtp.data(), rtp.size());

        Eigen::SparseMatrix<double, Eigen::RowMajor> A(dataSize, coeffs->size());
        A.setFromTriplets(TrL.begin(), TrL.end());

        Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
        SQR.compute(A);

        const auto X = SQR.solve(B);

        for (i = 0; i < coeffs->size(); i++) {
            (*coeffs)[i] = X(i);
        }

        ERR_RETURN
    }

    Err calcSpline(
            const QVector<double> &coeff,
            const QVector<double> &points,
            double x,
            double *y
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(coeff); ree;
        e = ErrorUtils::isNotEmpty(points); ree;

        double r = 0.0;
        int m = points.size();

        if (x <= points.at(0)) {
            r += coeff.at(0);
            r += coeff.at(1) * (x - points.at(0));
            *y = r;
            ERR_RETURN
        }
        else if (x >= points.at(m - 1)) {
            r += coeff.at((m - 1) * 2);
            r += coeff.at((m - 1) * 2 + 1) * (x - points.at(m - 1));
            *y = r;
        }

        int low = 0;
        int k = m - 1;
        while (k > low + 1) {
            int middle = (low + k) / 2;
            if (x < points.at(middle)) {
                k = middle;
            }
            else {
                low = middle;
            }
        }

        const double w = points.at(k) - points.at(k - 1);
        e = ErrorUtils::isFalse(MathUtils::tZero(w)); ree;

        const double u = (x - points.at(k - 1)) / w;
        const double v = u - 1.0;
        r += coeff.at((k - 1) * 2) * (1.0 + 2.0 * u) * std::pow(v, 2);
        r += coeff.at((k - 1) * 2 + 1) * w * u * std::pow(v, 2);

        r += coeff.at(k * 2) * std::pow(u, 2) * (1.0 - 2.0 * v);
        r += coeff.at(k * 2 + 1) * w * std::pow(u, 2) * v;

        *y = r;

        ERR_RETURN
    }


}//namespace
Err ScanTimeFromIRtMapper::mapRT(const QVector<QPair<double, double>> &data) {

    ERR_INIT

    const int minDataPoints = 3;
    e = ErrorUtils::isAboveThreshold(
            data.size(),
            minDataPoints,
            ErrorUtilsParam::IncludeThreshold
            ); ree;

    e = spline(data, m_segments, &m_coeffs, &m_points); ree;

    if (data.size() < minDataPoints) {
        ERR_RETURN
    }

    QVector<double> rtDiff;
    QVector<double> rtDiffSorted;
    QVector<QPair<double, double>> tempData;

    rtDiff.reserve(data.size());
    tempData.reserve(data.size());

    for (const QPair<double, double> &p : data) {

        double y;
        e = calcSpline(m_coeffs, m_points, p.first, &y); ree;
        rtDiff.push_back(std::abs(p.second - y));
    }

    rtDiffSorted = rtDiff;
    std::sort(rtDiffSorted.begin(), rtDiffSorted.end());

    const double max = rtDiffSorted[static_cast<int>(0.8 * rtDiffSorted.size())];
    for (int i = 0; i < data.size(); i++) {
        if (rtDiff[i] <= max + std::numeric_limits<double>::min()) {
            tempData.push_back(data[i]);
        }
    }

    spline(tempData, m_segments, &m_coeffs, &m_points);

    ERR_RETURN
}

Err ScanTimeFromIRtMapper::predictScanTime(double iRT, double *scanTime) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_coeffs); ree;
    e = calcSpline(m_coeffs, m_points, iRT, scanTime); ree;

    ERR_RETURN
}
