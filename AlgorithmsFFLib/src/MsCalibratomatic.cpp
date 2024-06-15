//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParquetReader.h"
#include "SqlUtils.h"
#include "XYMappermatic.h"

//#define USE_SPLINES

MsCalibratomatic::MsCalibratomatic()
: m_mzStDevMS1(-1.0)
, m_mzStDevMS2(-1.0)
, m_scanTimeStd(-1.0)
, m_isInitRT(false)
, m_isInitMS1(false)
, m_isInitMS2(false)
, m_polynomialOrderMassCal(3)
{}

namespace {

    Err generateMetricsIRTtoScanTime(
            const QVector<QPair<XVal, YVal>> &data,
            double *stDevScanTimeDiff
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(data); ree;

        const double testFraction = 0.2;
        const int seed = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;

        QVector<QPair<XVal, YVal>> trainingData;
        QVector<QPair<XVal, YVal>> testData;
        e = MathUtils::trainTestSplit(
                data,
                testFraction,
                seed,
                &trainingData,
                &testData
        ); ree;

        XYMappermatic mapperMetrics;
        e = mapperMetrics.init(trainingData); ree;

        QVector<QPair<double, double>> actualVsPredicted;
        QVector<double> diffs;
        for (const QPair<XVal, YVal> &pr : testData) {

            double yPredicted;
            e = mapperMetrics.predictY(pr.first, &yPredicted);
            actualVsPredicted.push_back({pr.second, yPredicted});
            diffs.push_back({pr.second - yPredicted});
        }

        const double rmse = MathUtils::rmse(actualVsPredicted);
        const double mean = MathUtils::mean(diffs);
        const double stDev = MathUtils::stDev(diffs);

        qDebug() << "iRT Cal Metrics: rmse" << rmse << "mean" << mean << "stDev" << stDev;

        *stDevScanTimeDiff = stDev;

        ERR_RETURN
    }

} //namespace
Err MsCalibratomatic::buildRTMapper(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;
    qDebug() << "Calibarating iRT->ScanTime";

    const auto insertLogic = [](const MsCalibarationReaderRow &r){
        return QPair(r.iRTPredicted, r.scanTime);
    };

    QVector<QPair<XVal, YVal>> data;
    data.reserve(msCalibarationReaderRows.size());
    std::transform(
            msCalibarationReaderRows.begin(),
            msCalibarationReaderRows.end(),
            std::back_inserter(data),
            insertLogic
            );

    double stDevScanTimeDiff;
    e = generateMetricsIRTtoScanTime(data, &stDevScanTimeDiff); ree;
    m_scanTimeStd = stDevScanTimeDiff;

    e = m_iRTtoScanTimeMapper.init(data); ree;
    e = ErrorUtils::isTrue(m_scanTimeStd > 0.0); ree;

    qDebug() << "scanTimeStDev" << m_scanTimeStd;
    m_isInitRT = true;

    ERR_RETURN
}

namespace {

    struct Inp {
        double mzSearched = -1.0;
        double mzFound = -1.0;
        double stDev = -1.0;
        double ppm = -1.0;
    };

    Err buildMzCalData(
            const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows,
            QVector<Inp> *inputs
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;

        int c = 0;
        for (const MsCalibarationReaderRow &row : msCalibarationReaderRows){

            const QVector<float> &mzSearchedVec = row.mzSearchedVec;
            const QVector<float> &mzFoundMeanVec = row.mzFoundMeanVec;
            const QVector<float> &mzFoundStDevVec = row.mzFoundStDevVec;

            e = ErrorUtils::isEqual(mzSearchedVec.size(), mzFoundMeanVec.size()); ree;
            e = ErrorUtils::isEqual(mzSearchedVec.size(), mzFoundStDevVec.size()); ree;

            for (int i = 0; i < mzFoundMeanVec.size(); i++) {
                Inp inp;
                inp.mzFound = mzFoundMeanVec.at(i);
                inp.mzSearched = mzSearchedVec.at(i);
                inp.stDev = std::max(mzFoundStDevVec.at(i), std::numeric_limits<float>::min());

                if (MathUtils::tZero(inp.mzSearched) || MathUtils::tZero(inp.mzFound)) {
                    continue;
                }

                inp.ppm = 1e6 * (inp.mzFound - inp.mzSearched) / inp.mzSearched;
                inputs->push_back(inp);
            }
        }

        qDebug() << inputs->size() << "points for mz recalibration";
        ERR_RETURN
    }

    Err removeMassOutliers(QVector<Inp> *inp) {

        ERR_INIT

        qDebug() << "Removing mz outliers for calibration";

        e = ErrorUtils::isNotEmpty(*inp); ree;

        QVector<double> ppms;
        ppms.reserve(inp->size());
        std::transform(
                inp->begin(),
                inp->end(),
                std::back_inserter(ppms),
                [](const Inp &i){return i.ppm;}
                );

        QVector<double> stDevs;
        stDevs.reserve(inp->size());
        std::transform(
                inp->begin(),
                inp->end(),
                std::back_inserter(stDevs),
                [](const Inp &i){return i.stDev;}
        );

        const double ppmMean = MathUtils::mean(ppms);
        const double ppmStDev = MathUtils::stDev(ppms);
        const auto ppmsMinMax = std::minmax_element(ppms.begin(), ppms.end());

        const double stDevMean = MathUtils::mean(stDevs);
        const double stDevMedian = MathUtils::median(stDevs);
        const double stDevStDev = MathUtils::stDev(stDevs);
        const auto stDevsMinMax = std::minmax_element(stDevs.begin(), stDevs.end());

        qDebug() << "ppmMean" << ppmMean << "ppmStDev"
                << ppmStDev << "ppmMin" << *ppmsMinMax.first << "ppmMax" << *ppmsMinMax.second;
        qDebug() << "stDevMean" << stDevMean << "stDevMean" << stDevStDev
                << "stDevMin" << *stDevsMinMax.first << "stDevMax" << *stDevsMinMax.second;

        constexpr int ppmStDevMultiplier = 3;
        const auto terminatorLogic = [&](const Inp &i){

            const double ppmTol = ppmStDev * ppmStDevMultiplier;
            const double ppmLo = ppmMean - ppmTol;
            const double ppmHi = ppmMean + ppmTol;
            const bool ppmOutOfRange = i.ppm < ppmLo || i.ppm > ppmHi;

            const bool stDevOutOfRange = i.stDev > stDevMedian;

            return ppmOutOfRange || stDevOutOfRange;
        };

        const auto terminator = std::remove_if(inp->begin(), inp->end(), terminatorLogic);

        inp->erase(terminator, inp->end());

        ERR_RETURN
    }

    Err generateMetricsMzReCal(
            const QVector<Inp> &inputs,
            double *stDevMz
    ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(inputs); ree;

        const double testFraction = 0.2;
        const int seed = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;

        QVector<QPair<double, double>> data;
        std::transform(
                inputs.begin(),
                inputs.end(),
                std::back_inserter(data),
                [](const Inp inp){return QPair<double, double>(inp.mzFound, inp.mzSearched);}
        );

        QVector<QPair<XVal, YVal>> trainingData;
        QVector<QPair<XVal, YVal>> testData;
        e = MathUtils::trainTestSplit(
                data,
                testFraction,
                seed,
                &trainingData,
                &testData
        ); ree;

        QVector<double> mzToRecalCalCurveCoeffs;
        XYMappermatic mapperMetrics;

        Eigen::MatrixX<double> mat(trainingData.size(), 2);
        for (int i = 0; i < trainingData.size(); i++) {
            const QPair<XVal, YVal> &pr = trainingData.at(i);
            mat.coeffRef(i, 0) = pr.first;
            mat.coeffRef(i, 1) = pr.second;
        }
        const int polynomialOrder = 5;
        EigenUtils::fitPolynomialQRDecomposition(mat, polynomialOrder, &mzToRecalCalCurveCoeffs);

        QVector<QPair<double, double>> actualVsPredicted;

        QVector<double> ppmOriginals;
        QVector<double> ppmReCals;
        for (const QPair<XVal, YVal> &pr : testData) {

            double yPredicted = 0.0;

            for (int i = 0; i < mzToRecalCalCurveCoeffs.size(); i++) {
                yPredicted += mzToRecalCalCurveCoeffs.at(i) * std::pow(pr.first, i);
            }

            actualVsPredicted.push_back({pr.second, yPredicted});

            const double ppmOriginal = 1e6 * (pr.first - pr.second) / pr.second;
            const double ppmReCal = 1e6 * (yPredicted - pr.second) / pr.second;

            ppmOriginals.push_back(ppmOriginal);
            ppmReCals.push_back(ppmReCal);

        }

        const double rmse = MathUtils::rmse(actualVsPredicted);
        const double meanOriginal = MathUtils::mean(ppmOriginals);
        const double stDevOriginal = MathUtils::stDev(ppmOriginals);
        const double meanReCal = MathUtils::mean(ppmReCals);
        const double stDevReCal = MathUtils::stDev(ppmReCals);

        qDebug() << "Mz Cal Metrics: rmse" << rmse;
        qDebug() << "meanOG PPM" << meanOriginal << "stDevOG PPM" << stDevOriginal;
        qDebug() << "meanReCal PPM" << meanReCal << "stDevReCal PPM" << stDevReCal;

        *stDevMz = stDevReCal;

        ERR_RETURN
    }

    Err buildMzCalibrationCurve(
        const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows,
        const int polynomialOrder,
        double *mzStDev,
        QVector<double> *calibrationCurveCoEffs
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;

        //TODO figure out how to tighten this up.  StDev doesn't change implying that this algo is simply shifting everything.

        QVector<Inp> inputs;
        e = buildMzCalData(msCalibarationReaderRows, &inputs); ree;
        e = removeMassOutliers(&inputs); ree;

        e = generateMetricsMzReCal(inputs, mzStDev); ree;

        Eigen::MatrixX<double> mat(inputs.size(), 2);
        for (int i = 0; i < inputs.size(); i++) {
            const Inp &inp = inputs.at(i);
            mat.coeffRef(i, 0) = inp.mzFound;
            mat.coeffRef(i, 1) = inp.mzSearched;
        }
        EigenUtils::fitPolynomialQRDecomposition(mat, polynomialOrder, calibrationCurveCoEffs);

        ERR_RETURN
    }

}//namespace
Err MsCalibratomatic::setMassCalibrationCoeffs(
    const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows,
    const MSLevelEnum msLevelEnum
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;

    if (msLevelEnum == MSLevelEnum::MS1) {

        e = buildMzCalibrationCurve(
            msCalibarationReaderRows,
            m_polynomialOrderMassCal,
            &m_mzStDevMS1,
            &m_calibrationCurveCoeffsMS1
        ); ree;

        qDebug() << "mzStDevMS1" << m_mzStDevMS1;
        e = ErrorUtils::isTrue(m_mzStDevMS1 > 0.0); ree;

        m_isInitMS1 = true;

        ERR_RETURN
    }

    e = buildMzCalibrationCurve(
        msCalibarationReaderRows,
        m_polynomialOrderMassCal,
        &m_mzStDevMS2,
        &m_calibrationCurveCoeffsMS2
    ); ree;

    qDebug() << "mzStDevMS2" << m_mzStDevMS2;
    e = ErrorUtils::isTrue(m_mzStDevMS2 > 0.0); ree;

    m_isInitMS2 = true;

    ERR_RETURN
}

Err MsCalibratomatic::recalibrateScanPoints(
        const MSLevelEnum &msLevel,
        const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
        ) {

    ERR_INIT

    if (msLevel == MSLevelEnum::MS2) {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS2); ree;
    }
    else {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS1); ree;
    }

    const QVector<double> calibrationCoeffs = msLevel == MSLevelEnum::MS2 ? m_calibrationCurveCoeffsMS2 : m_calibrationCurveCoeffsMS1;

    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); it++) {

        ScanPoints *scanPoints = it.value();

        for (ScanPoint &sp : *scanPoints) {

            double mzRecal = 0.0;
            for (int i = 0; i < calibrationCoeffs.size(); i++) {
                mzRecal += calibrationCoeffs.at(i) * std::pow(sp.x(), i);
            }
        }
    }

    ERR_RETURN
}

Err MsCalibratomatic::recalibrateScanPoints(
        const MSLevelEnum &msLevel,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
){

    ERR_INIT

    if (msLevel == MSLevelEnum::MS2) {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS2); ree;
    }
    else {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS1); ree;
    }

    const QVector<double> calibrationCoeffs = msLevel == MSLevelEnum::MS2 ? m_calibrationCurveCoeffsMS2 : m_calibrationCurveCoeffsMS1;

    for (auto it = scanNumberVsScanPoints->begin(); it != scanNumberVsScanPoints->end(); it++) {

        ScanPoints &scanPoints = it.value();

        for (ScanPoint &sp : scanPoints) {
            double mzRecal = 0.0;

            for (int i = 0; i < calibrationCoeffs.size(); i++) {
                mzRecal += calibrationCoeffs.at(i) * std::pow(sp.x(), i);
            }

            sp.rx() = static_cast<float>(mzRecal);

        }
    }

    ERR_RETURN
}

Err MsCalibratomatic::recalibrateScanPoint(
    const MSLevelEnum& msLevel,
    float mzVal,
    float* mzValRecal
    ) {

    ERR_INIT

    if (msLevel == MSLevelEnum::MS2) {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS2); ree;
    }
    else {
        e = ErrorUtils::isNotEmpty(m_calibrationCurveCoeffsMS1); ree;
    }

    const QVector<double> calibrationCoeffs = msLevel == MSLevelEnum::MS2 ? m_calibrationCurveCoeffsMS2 : m_calibrationCurveCoeffsMS1;

    double mzRecal = 0.0f;
    for (int i = 0; i < calibrationCoeffs.size(); i++) {
        mzRecal += calibrationCoeffs.at(i) * std::pow(mzVal, i);
    }
    *mzValRecal = static_cast<float>(mzRecal);

    ERR_RETURN
}

float MsCalibratomatic::mzStDevMS1() const {
    return m_mzStDevMS1;
}

float MsCalibratomatic::mzStDevMS2() const {
    return m_mzStDevMS2;
}

float MsCalibratomatic::scanTimeStDev(float nStdDevs /* = 1.0 */) const {
    return static_cast<float>(m_scanTimeStd) * nStdDevs;
}

Err MsCalibratomatic::predictScanTime(float iRT, float *predictedScanTime) const {
    ERR_INIT

    double predictedScanTimeDouble;
    e = m_iRTtoScanTimeMapper.predictY(static_cast<double>(iRT), &predictedScanTimeDouble); ree;
    *predictedScanTime = static_cast<float>(predictedScanTimeDouble);

    ERR_RETURN
}

bool MsCalibratomatic::isInitRT() const {
    return m_isInitRT;
}

bool MsCalibratomatic::isInitCalMS1() const {
    return m_isInitMS1;
}

bool MsCalibratomatic::isInitCalMS2() const {
    return m_isInitMS2;
}

void MsCalibratomatic::setScanTimeStDev(double val) {
    m_scanTimeStd = val;
    qDebug() << "scanTimeStDev has been set to:" << scanTimeStDev() << "seconds";
    qDebug() << "scanTimeStDev x 3:" << scanTimeStDev(3) << "seconds";
}

void MsCalibratomatic::setMzStDevMS2(double val) {
    m_mzStDevMS2 = val;
    qDebug() << "MzStdDevMS2 has been set to:" << m_mzStDevMS2;
}
