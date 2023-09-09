//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "ErrorUtils.h"
#include "MsCalibrationReader.h"
#include "ParquetReader.h"
#include "SqlUtils.h"
#include "XYMappermatic.h"


MsCalibratomatic::MsCalibratomatic()
: m_stDevNew(-1.0)
, m_scanTimeWindowNew(-1.0)
{}

Err MsCalibratomatic::init(
        const PythiaParameters &pythiaParameters,
        const QString &msCalibrationFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree
    m_params = pythiaParameters;

    e = ErrorUtils::fileExists(msCalibrationFilePath); ree
    m_msCalibrationFilePath = msCalibrationFilePath;

    e = buildIRTCalibrator(); ree
    e = buildMzCalibrator(); ree

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
        const int seed = 666;

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

        *stDevScanTimeDiff = stDev * S_GLOBAL_SETTINGS.STDEV_MULTIPLIER;

        ERR_RETURN
    }

} //namespace
Err MsCalibratomatic::buildIRTCalibrator() {

    ERR_INIT

    qDebug() << "Calibarating iRT->ScanTime";

    QVector<MsCalibarationReaderRow> msCalibarationReaderRows;
    e = ParquetReader::read(
            m_msCalibrationFilePath,
            &msCalibarationReaderRows
            ); ree;

    const auto insertLogic = [](const MsCalibarationReaderRow &r){
        return QPair(r.iRTPredicted, r.scanTime);
    };

    QVector<QPair<XVal, YVal>> data;
    std::transform(
            msCalibarationReaderRows.begin(),
            msCalibarationReaderRows.end(),
            std::back_inserter(data),
            insertLogic
            );

    double stDevScanTimeDiff;
    e = generateMetricsIRTtoScanTime(data, &stDevScanTimeDiff); ree;
    e = setScanTimeWindowNew(stDevScanTimeDiff); ree;

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

            const QVector<double> &mzSearchedVec = row.mzSearchedVec;
            const QVector<double> &mzFoundMeanVec = row.mzFoundMeanVec;
            const QVector<double> &mzFoundStDevVec = row.mzFoundStDevVec;

            e = ErrorUtils::isEqual(mzSearchedVec.size(), mzFoundMeanVec.size()); ree;
            e = ErrorUtils::isEqual(mzSearchedVec.size(), mzFoundStDevVec.size()); ree;

            for (int i = 0; i < mzFoundMeanVec.size(); i++) {
                Inp inp;
                inp.mzFound = mzFoundMeanVec.at(i);
                inp.mzSearched = mzSearchedVec.at(i);
                inp.stDev = std::max(mzFoundStDevVec.at(i), std::numeric_limits<double>::min());

                if (MathUtils::tZero(inp.mzSearched) || MathUtils::tZero(inp.mzFound)) {
                    continue;
                }

                inp.ppm = 1e6 * (inp.mzFound - inp.mzSearched) / inp.mzSearched;
                inputs->push_back(inp);
            }
        }

        ERR_RETURN
    }

    Err removeMassOutliers(QVector<Inp> *inp) {

        ERR_INIT

        qDebug() << "Removing mz outliers for calibration";

        e = ErrorUtils::isNotEmpty(*inp); ree;

        QVector<double> ppms;
        std::transform(
                inp->begin(),
                inp->end(),
                std::back_inserter(ppms),
                [](const Inp &i){return i.ppm;}
                );

        QVector<double> stDevs;
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
        const int seed = 666;

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

        *stDevMz = stDevReCal * S_GLOBAL_SETTINGS.STDEV_MULTIPLIER;

        ERR_RETURN
    }

}//namespace
Err MsCalibratomatic::buildMzCalibrator() {

    ERR_INIT

    //TODO figure out how to tighten this up.  StDev doesn't change implying that this algo is simply shifting everything.

    QVector<MsCalibarationReaderRow> msCalibarationReaderRows;
    e = ParquetReader::read(
            m_msCalibrationFilePath,
            &msCalibarationReaderRows
    ); ree;

    QVector<Inp> inputs;
    e = buildMzCalData(msCalibarationReaderRows, &inputs); ree;
    e = removeMassOutliers(&inputs); ree;

    double stDevMz;
    e = generateMetricsMzReCal(inputs, &stDevMz); ree;

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
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        QMap<ScanNumber, ScanPoints> *reCalScanNumberVsScanPoints
        ) {

    ERR_INIT

    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPoints = it.value();

        ScanPoints scanPointsRecal;
        for (const ScanPoint &sp : scanPoints) {
            double mzRecal;
            e = m_mzToRecalMz.predictY(sp.x(), &mzRecal); ree;
            scanPointsRecal.push_back({mzRecal, sp.y()});
        }

        reCalScanNumberVsScanPoints->insert(scanNumber, scanPointsRecal);
    }

    ERR_RETURN
}

double MsCalibratomatic::newStDev() {
    return m_stDevNew;
}

Err MsCalibratomatic::setScanTimeWindowNew(double scanTimeWindow) {

    ERR_INIT

    e = ErrorUtils::isFalse(MathUtils::tZero(scanTimeWindow)); ree;
    e = ErrorUtils::isTrue(scanTimeWindow > 0.0); ree;

    m_scanTimeWindowNew = scanTimeWindow;
    qDebug() << "scanTimeWindowNew set to" << m_scanTimeWindowNew;

    ERR_RETURN
}

Err MsCalibratomatic::predictScanTime(double iRT, double *predictedScanTime) {
    ERR_INIT
    e = m_iRTtoScanTimeMapper.predictY(iRT, predictedScanTime); ree;
    ERR_RETURN
}
