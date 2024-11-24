//
// Created by anichols on 11/07/2021.
//

#include "MsReaderMzMLLazyLoad.h"
#include "FastaReader.h"
#include "ParallelUtils.h"

#include <QString>
#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>


class MsReaderMZMLLazyLoadTests : public QObject
{
    Q_OBJECT

public:
    MsReaderMZMLLazyLoadTests() = default;
    ~MsReaderMZMLLazyLoadTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    static void openFileTest();
    static void troubleShoot();

};


void MsReaderMZMLLazyLoadTests::openFileTest() {

    QSKIP("Turn on after complete and making tests.");

    // const QString &msParquetFilePath = QDir(qApp->applicationDirPath()).filePath("1min.mzML");
    //
    // ERR_INIT
    //
    // QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));
    //
    // MsReaderMzML reader;
    // e = reader.openFile(msParquetFilePath);
    // QCOMPARE(e, Error::eNoError);
    // QCOMPARE(reader.m_msScanInfo.size(), 372);
    // QCOMPARE(reader.m_scanPoints.size(), 372);
    //
    // MsScanInfo msScanInfo;
    // e = reader.getMsScanInfo(372, &msScanInfo);
    // QCOMPARE(e, eNoError);
    // QCOMPARE(msScanInfo.scanNumber, 372);
    // QCOMPARE(msScanInfo.msLevel, 1);
    // QCOMPARE(msScanInfo.collisionEnergy, -1);
    // QCOMPARE(QString::number(msScanInfo.scanTime), "10.0411");
    // QCOMPARE(QString::number(msScanInfo.precursorTargetMz), "-1");
    // QCOMPARE(msScanInfo.isoWindowLower, -1);
    // QCOMPARE(msScanInfo.isoWindowUpper, -1);
    //
    // QPair<Err, ScanPoints*> scanPointsResult = reader.getScanPoints(msScanInfo.scanNumber);
    // e = scanPointsResult.first;
    // QCOMPARE(e, eNoError);
    //
    // const ScanPoints scanPoints = *scanPointsResult.second;
    //
    // QCOMPARE(int(scanPoints.first().x()), 150);
    // QCOMPARE(int(scanPoints.first().y()), 814);
    // QCOMPARE(int(scanPoints.last().x()), 1250);
    // QCOMPARE(int(scanPoints.last().y()), 2646);
    // QCOMPARE(scanPoints.size(), 2092);

}

namespace {
    Err parallelLogic(
        const MsScanInfo &msi,
        MsReaderMzMLLazyLoad *msReaderMzMLLazyLoad
        ) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isTrue(msReaderMzMLLazyLoad->isInit()); ree;

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
        e = msReaderMzMLLazyLoad->getMzTargetScanPoints(msi.targetKey(), &scanNumberVsScanPoints); ree;
        qDebug() << msi.targetKey() << scanNumberVsScanPoints.size() << "scanNumberVsScanPoints Size yo!" << et.elapsed();
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "ms";

        ERR_RETURN
    }

    QVector<Err> parallelLogicLoop(
        const QVector<MsScanInfo> &msis,
        const QString &filePath
        ) {

        ERR_INIT

        QVector<Err> retVal;

        MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
        e = msReaderMzMLLazyLoad.openFile(filePath);

        for (const MsScanInfo &msi : msis) {
            e = parallelLogic(msi, &msReaderMzMLLazyLoad);
            retVal.push_back(e);
        }

        return retVal;
    }

}//namespace
void MsReaderMZMLLazyLoadTests::troubleShoot() {

    ERR_INIT

    QSKIP("Troubleshooting");
    //TODO make more tests.

    const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
    // const QString filename = "/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";
    // const QString filename = "/home/anichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML";

    // const QString &filename = QDir(qApp->applicationDirPath()).filePath("1min.mzML");

    MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
    e = msReaderMzMLLazyLoad.openFile(filename);
    QCOMPARE(e, eNoError);
    QCOMPARE(msReaderMzMLLazyLoad.m_mzTargetVsScanInfosPntrs.size(), 201);

    QVector<MsScanInfo> scanInfos = msReaderMzMLLazyLoad.getUniqueTandemMsScanInfos();

#define RUN_PARALLEL_TEST
#ifdef RUN_PARALLEL_TEST

    QVector<QVector<MsScanInfo>> msScanInfosesTranced;
    e = ParallelUtils::trancheVectorForParallelization(
        scanInfos,
        ParallelUtils::numberOfAvailableSystemProcessors(),
        &msScanInfosesTranced
        );
    QCOMPARE(e, eNoError);

    const auto trainingLogicBinder = std::bind(
        parallelLogic,
        std::placeholders::_1,
        &msReaderMzMLLazyLoad
        );

    QFuture<Err> futures = QtConcurrent::mapped(
        scanInfos,
        trainingLogicBinder
        );
    futures.waitForFinished();

#else
    for (const MsScanInfo& msi : scanInfos) {
        e = parallelLogic(msi, &msReaderMzMLLazyLoad);
        QCOMPARE(e, eNoError);
    }
#endif
}


QTEST_MAIN(MsReaderMZMLLazyLoadTests)
#include "MsReaderMZMLLazyLoadTests.moc"
