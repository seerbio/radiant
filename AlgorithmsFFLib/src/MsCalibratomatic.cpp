//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParquetReader.h"
#include "SqlUtils.h"
#include "XYMappermatic.h"


MsCalibratomatic::MsCalibratomatic()
: m_mzStDev(-1.0)
, m_scanTimeStd(-1.0)
, m_isInit(false)
{}

Err MsCalibratomatic::init(const QString &msCalibrationFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msCalibrationFilePath); ree
    m_msCalibrationFilePath = msCalibrationFilePath;

    e = ParquetReader::read(
            m_msCalibrationFilePath,
            &m_msCalibarationReaderRows
    ); ree;

    e = init(m_msCalibarationReaderRows); ree;

    ERR_RETURN
}


Err MsCalibratomatic::init(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;
    m_msCalibarationReaderRows = msCalibarationReaderRows;

    e = buildIRTCalibrator(); ree
    e = buildMzCalibrator(); ree

    e = ErrorUtils::isTrue(m_mzStDev > 0.0); ree;
    e = ErrorUtils::isTrue(m_scanTimeStd > 0.0); ree;

    qDebug() << "mzStDev" << m_mzStDev;
    qDebug() << "scanTimeStDev" << m_scanTimeStd;

    m_isInit = true;

    ERR_RETURN
}

Err MsCalibratomatic::initRtOnly(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows) {
    ERR_INIT

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;
    m_msCalibarationReaderRows = msCalibarationReaderRows;

    e = buildIRTCalibrator(); ree
    e = ErrorUtils::isTrue(m_scanTimeStd > 0.0); ree;

    qDebug() << "scanTimeStDev" << m_scanTimeStd;
    m_isInit = true;

    ERR_RETURN
}

Err MsCalibratomatic::initMzOnly(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msCalibarationReaderRows); ree;
    m_msCalibarationReaderRows = msCalibarationReaderRows;

    e = buildMzCalibrator(); ree

    e = ErrorUtils::isTrue(m_mzStDev > 0.0); ree;

    qDebug() << "mzStDev" << m_mzStDev;
    m_isInit = true;

    ERR_RETURN
}

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
Err MsCalibratomatic::buildIRTCalibrator() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msCalibarationReaderRows); ree;
    qDebug() << "Calibarating iRT->ScanTime";

    const auto insertLogic = [](const MsCalibarationReaderRow &r){
        return QPair(r.iRTPredicted, r.scanTime);
    };

    QVector<QPair<XVal, YVal>> data;
    data.reserve(m_msCalibarationReaderRows.size());
    std::transform(
            m_msCalibarationReaderRows.begin(),
            m_msCalibarationReaderRows.end(),
            std::back_inserter(data),
            insertLogic
            );

    double stDevScanTimeDiff;
    e = generateMetricsIRTtoScanTime(data, &stDevScanTimeDiff); ree;
    m_scanTimeStd = stDevScanTimeDiff;

    e = m_iRTtoScanTimeMapper.init(data); ree;

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
        const double stDevStDev = MathUtils::stDev(stDevs);
        const auto stDevsMinMax = std::minmax_element(stDevs.begin(), stDevs.end());

        qDebug() << "ppmMean" << ppmMean << "ppmStDev"
                << ppmStDev << "ppmMin" << *ppmsMinMax.first << "ppmMax" << *ppmsMinMax.second;
        qDebug() << "stDevMean" << stDevMean << "stDevMean" << stDevStDev
                << "stDevMin" << *stDevsMinMax.first << "stDevMax" << *stDevsMinMax.second;

        const auto terminatorLogic = [&](const Inp &i){

            const double ppmTol = ppmStDev * S_GLOBAL_SETTINGS.STDEV_MULTIPLIER;
            const double ppmLo = ppmMean - ppmTol;
            const double ppmHi = ppmMean + ppmTol;
            const bool ppmOutOfRange = i.ppm < ppmLo || i.ppm > ppmHi;

            const double stDevTol = stDevStDev * (S_GLOBAL_SETTINGS.STDEV_MULTIPLIER);
            const double stDevLo = stDevMean - stDevTol;
            const double stDevHi = stDevMean + stDevTol;
            const bool stDevOutOfRange = i.stDev < stDevLo || i.stDev > stDevHi;

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

        XYMappermatic mapperMetrics;
        e = mapperMetrics.init(trainingData); ree;

        QVector<QPair<double, double>> actualVsPredicted;

        QVector<double> ppmOriginals;
        QVector<double> ppmReCals;
        for (const QPair<XVal, YVal> &pr : testData) {

            double yPredicted;
            e = mapperMetrics.predictY(pr.first, &yPredicted);
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

}//namespace
Err MsCalibratomatic::buildMzCalibrator() {

    ERR_INIT


    e = ErrorUtils::isNotEmpty(m_msCalibarationReaderRows); ree;

    //TODO figure out how to tighten this up.  StDev doesn't change implying that this algo is simply shifting everything.

    QVector<Inp> inputs;
    e = buildMzCalData(m_msCalibarationReaderRows, &inputs); ree;
    e = removeMassOutliers(&inputs); ree;

    double stDevMz;
    e = generateMetricsMzReCal(inputs, &stDevMz); ree;
    m_mzStDev = stDevMz;

    QVector<QPair<double, double>> data;
    std::transform(
            inputs.begin(),
            inputs.end(),
            std::back_inserter(data),
            [](const Inp inp){return QPair<double, double>(inp.mzFound, inp.mzSearched);}
    );

    e = m_mzToRecalMz.init(data); ree;

    ERR_RETURN
}

Err MsCalibratomatic::recalibrateScanPoints(
        const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
        ) {

    ERR_INIT

    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); it++) {

        ScanPoints *scanPoints = it.value();

        for (ScanPoint &sp : *scanPoints) {
            double mzRecal;
            e = m_mzToRecalMz.predictY(sp.x(), &mzRecal); ree;
            sp.rx() = mzRecal;
        }
    }

    ERR_RETURN
}


Err MsCalibratomatic::recalibrateScanPoints(
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
) {

    ERR_INIT

    for (auto it = scanNumberVsScanPoints->begin(); it != scanNumberVsScanPoints->end(); it++) {

        ScanPoints &scanPoints = it.value();

        for (ScanPoint &sp : scanPoints) {
            double mzRecal;
            e = m_mzToRecalMz.predictY(sp.x(), &mzRecal); ree;
            sp.rx() = mzRecal;
        }
    }

    ERR_RETURN
}

float MsCalibratomatic::mzStDev() {
    return m_mzStDev;
}

float MsCalibratomatic::scanTimeStDev(int nStdDevs /* = 1 */) {
    return m_scanTimeStd * nStdDevs;
}

Err MsCalibratomatic::predictScanTime(float iRT, float *predictedScanTime) const {
    ERR_INIT

    double predictedScanTimeDouble;
    e = m_iRTtoScanTimeMapper.predictY(static_cast<double>(iRT), &predictedScanTimeDouble); ree;
    *predictedScanTime = static_cast<float>(predictedScanTimeDouble);
    ERR_RETURN
}

bool MsCalibratomatic::isInit() const {
    return m_isInit;
}

void MsCalibratomatic::setScanTimeStDev(double val) {
    m_scanTimeStd = val;
    qDebug() << "scanTimeStDev has been set to:" << scanTimeStDev() << "seconds";
    qDebug() << "scanTimeStDev x 3:" << scanTimeStDev(3) << "seconds";
}
