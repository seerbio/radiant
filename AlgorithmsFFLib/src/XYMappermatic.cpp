//
// Created by anichols on 8/24/23.
//

#include "XYMappermatic.h"

#include "ErrorUtils.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <iostream>

XYMappermatic::XYMappermatic()
: m_segments(0)
, m_xSegments(20)
, m_minXPredBin(20)
{}


Err XYMappermatic::init(const QString &iRTRecalibrationFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(iRTRecalibrationFilePath); ree;

    QVector<XYMappermaticRow> xyMappermaticRows;
    e  = CSVReader::read(
            iRTRecalibrationFilePath,
            S_GLOBAL_SETTINGS.COMMA,
            &xyMappermaticRows
    ); ree;

    qDebug() << xyMappermaticRows.size() << "X Rows";

    const auto insertLogic
        = [](const XYMappermaticRow &r){return QPair<double, double>(r.x, r.y);};

    QVector<QPair<XVal, YVal>> data;
    data.reserve(xyMappermaticRows.size());
    std::transform(
            xyMappermaticRows.begin(),
            xyMappermaticRows.end(),
            std::back_inserter(data),
            insertLogic
            );

    e = init(data); ree;

    ERR_RETURN
}

Err XYMappermatic::init(const QVector<QPair<XVal, YVal>> &_data) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_data); ree;

    QVector<QPair<XVal, YVal>> data = _data;

    const auto sortLogicXAsc
            = [](const QPair<XVal , YVal> &l, const QPair<double, double> &r){return l.first < r.first;};

    std::sort(data.begin(), data.end(), sortLogicXAsc);

    m_segments = std::min(
            m_xSegments,
            static_cast<int>(std::max(1.0, 2.0 * sqrt(data.size() / static_cast<double>(m_minXPredBin))))
    );

    e = mapXtoY(data); ree;

    ERR_RETURN
}


namespace {

    Err spline(
            const QVector<QPair<XVal , YVal>> &data,
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

        rtp.resize(dataSize);
        rte.resize(dataSize);

        m = segments;

        coeffs->resize(2 * m);

        rts.resize(m + 1);
        points->resize(m);

        for (i = 0; i < dataSize; i++) {
            rtp[i] = data[i].second;
            rte[i] = data[i].first;
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
                    (*points)[k] = points->at(k + 1);
                }
                m--;
            }
        }

        TrL.clear();

        for (i = k = 0; i < dataSize; i++) {
            double x = rte[i];

            if (!k && x <= points->at(0)) {

                TrL.emplace_back(i, 0, 1.0);
                TrL.emplace_back(i, 1, static_cast<double>(x - points->at(0)));
            }
            else {
                for (; k < m && points->at(k) < x; k++);
                if (k == m) {
                    TrL.emplace_back(i, (k - 1) * 2, 1.0);
                    TrL.emplace_back(i, (k - 1) * 2 + 1, static_cast<double>(x - points->at(k - 1)));
                }
                else {

                    w = points->at(k) - points->at(k - 1);
                    e = ErrorUtils::isFalse(MathUtils::tZero(w)); ree

                    const double u = (x - points->at(k - 1)) / w;
                    const double v = u - 1.0;

                    TrL.emplace_back(i, (k - 1) * 2, static_cast<double>((1.0 + 2.0 * u) * std::pow(v, 2)));
                    TrL.emplace_back(i, (k - 1) * 2 + 1, static_cast<double>(w * u * std::pow(v, 2)));
                    TrL.emplace_back(i, k * 2, static_cast<double>(std::pow(u, 2) * (1.0 - 2.0 * v)));
                    TrL.emplace_back(i, k * 2 + 1, static_cast<double>(w * std::pow(u, 2) * v));
                }
            }
        }

        const Eigen::VectorXd B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rtp.data(), rtp.size());

        Eigen::SparseMatrix<double, Eigen::RowMajor> A(dataSize, coeffs->size());
        for (const Triplet &t : TrL) {

            const bool coorIsOutOfRange = t.row() >= dataSize || t.col() >= coeffs->size() || t.row() < 0 || t.col() < 0;
            if (coorIsOutOfRange) {
                qDebug() << t.row() << t.col() << dataSize << coeffs->size() << "OOR Error";
            }

            e = ErrorUtils::isFalse(coorIsOutOfRange); ree;
            A.coeffRef(t.row(), t.col()) = t.value();
        }

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
            ERR_RETURN
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
Err XYMappermatic::mapXtoY(const QVector<QPair<XVal, YVal>> &data) {

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

    QVector<double> yDiff;
    QVector<double> yDiffSorted;
    QVector<QPair<XVal, YVal>> tempData;

    yDiff.reserve(data.size());
    tempData.reserve(data.size());

    for (const QPair<XVal, YVal> &p : data) {
        double y;
        e = calcSpline(m_coeffs, m_points, p.first, &y); ree;
        yDiff.push_back(std::abs(p.second - y));
    }

    yDiffSorted = yDiff;
    std::sort(yDiffSorted.begin(), yDiffSorted.end());

    const double max = yDiffSorted[static_cast<int>(0.8 * yDiffSorted.size())];
    for (int i = 0; i < data.size(); i++) {
        if (yDiff[i] <= max + std::numeric_limits<double>::min()) {
            tempData.push_back(data[i]);
        }
    }

    spline(tempData, m_segments, &m_coeffs, &m_points);

    ERR_RETURN
}

Err XYMappermatic::predictY(double x, double *y) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_coeffs); ree;
    e = calcSpline(m_coeffs, m_points, x, y); ree;

    ERR_RETURN
}

Err XYMappermatic::_splineTestAccess(
        const QVector<QPair<XVal, YVal>> &data,
        int segments,
        QVector<double> *coeffs,
        QVector<double> *points
        ) {

    ERR_INIT
    e = spline(data, segments, coeffs, points); ree
    ERR_RETURN
}
