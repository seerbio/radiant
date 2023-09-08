//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "MsCalibrationReader.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "MsFrameScoretronProcessormatic.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "PSMsReader.h"
#include "TandemFragmentPredictotron.h"
#include "XYMappermatic.h"

#include <Eigen/Dense>

#include <QtConcurrent/QtConcurrent>


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

    Err generateMetrics(
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

        qDebug() << "Cal Metrics: rmse" << rmse << "mean" << mean << "stDev" << stDev;

        *stDevScanTimeDiff = stDev * S_GLOBAL_SETTINGS.STDEV_MULTIPLIER;

        ERR_RETURN
    }

} //namespace
Err MsCalibratomatic::buildIRTCalibrator() {

    ERR_INIT

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
    e = generateMetrics(data, &stDevScanTimeDiff); ree;
    e = setScanTimeWindowNew(stDevScanTimeDiff); ree;

    e = m_iRTtoScanTimeMapper.init(data); ree;

    ERR_RETURN
}

Err MsCalibratomatic::buildMzCalibrator() {

    ERR_INIT

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


    ERR_RETURN
}


namespace {


}//namespace
Err MsCalibratomatic::recalibrateScanPoints(
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        QMap<ScanNumber, ScanPoints> *recalScanNumberVsScanPoints
        ) {

    ERR_INIT



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
