#include <QtTest>
#include <QCoreApplication>

#include "DeMultiplexScanerator.h"
#include "MsReaderParquet.h"

#include <QElapsedTimer>

#include "Eigen/Dense"
#include <unsupported/Eigen/NNLS>

#include <iostream>

class DeMultiplexScaneratorTests : public QObject {
Q_OBJECT

public:
    DeMultiplexScaneratorTests() = default;

    ~DeMultiplexScaneratorTests() override = default;

private slots:

    static void buildTransitionsMatrixTest();

    static void buildScanMaskMatrixTest();

    static void deMultiplexScansTest();

};

void DeMultiplexScaneratorTests::buildTransitionsMatrixTest() {

    ERR_INIT

    const QString prqFFFile = "/home/anichols/Downloads/TESTING_20240130RC6_30min_iso-4-500-510-slide-7test-15k.mzML.prqFF";

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(prqFFFile);
    QCOMPARE(e, eNoError);

    const ScanNumber scanNumber = 21797;
    const QPair<Err, ScanPoints*> scan1 = msReaderParquet.getScanPoints(scanNumber);
    QCOMPARE(scan1.first, eNoError);

    const QPair<Err, ScanPoints*> scan2 = msReaderParquet.getScanPoints(scanNumber + 1);
    QCOMPARE(scan2.first, eNoError);

    const QPair<Err, ScanPoints*> scan3 = msReaderParquet.getScanPoints(scanNumber + 2);
    QCOMPARE(scan3.first, eNoError);

    const QPair<Err, ScanPoints*> scan4 = msReaderParquet.getScanPoints(scanNumber + 3);
    QCOMPARE(scan4.first, eNoError);

    const QVector<ScanPoints*> scans = {
            scan1.second,
            scan2.second,
            scan3.second,
            scan4.second
    };

    DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);

    QElapsedTimer et;
    et.start();
    e =  deMultiplexScanerator._buildTransitionMatrixTestAccess(scans);
    qDebug() << et.elapsed();
    QCOMPARE(e, eNoError);

}

void DeMultiplexScaneratorTests::buildScanMaskMatrixTest() {

    ERR_INIT

    const QString prqFFFile = "/home/anichols/Downloads/TESTING_20240130RC6_30min_iso-4-500-510-slide-7test-15k.mzML.prqFF";

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(prqFFFile);
    QCOMPARE(e, eNoError);

    const ScanNumber scanNumber = 21797;

    MsScanInfo msScanInfo1;
    e = msReaderParquet.getMsScanInfo(scanNumber, &msScanInfo1);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo2;
    e = msReaderParquet.getMsScanInfo(scanNumber + 1, &msScanInfo2);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo3;
    e = msReaderParquet.getMsScanInfo(scanNumber + 2, &msScanInfo3);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo4;
    e = msReaderParquet.getMsScanInfo(scanNumber + 3, &msScanInfo4);
    QCOMPARE(e, eNoError);

    const QVector<MsScanInfo> msScanInfos = {
            msScanInfo1,
            msScanInfo2,
            msScanInfo3,
            msScanInfo4
    };

    DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
    e = deMultiplexScanerator._buildScanMaskMatrixTestAccess(msScanInfos);
    QCOMPARE(e, eNoError);

}

void DeMultiplexScaneratorTests::deMultiplexScansTest() {

    ERR_INIT

    const QString prqFFFile = "/home/anichols/Downloads/TESTING_20240130RC6_30min_iso-4-500-510-slide-7test-15k.mzML.prqFF";

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(prqFFFile);
    QCOMPARE(e, eNoError);

    const ScanNumber scanNumber = 21797;
    const QPair<Err, ScanPoints*> scan1 = msReaderParquet.getScanPoints(scanNumber);
    QCOMPARE(scan1.first, eNoError);

    const QPair<Err, ScanPoints*> scan2 = msReaderParquet.getScanPoints(scanNumber + 1);
    QCOMPARE(scan2.first, eNoError);

    const QPair<Err, ScanPoints*> scan3 = msReaderParquet.getScanPoints(scanNumber + 2);
    QCOMPARE(scan3.first, eNoError);

    const QPair<Err, ScanPoints*> scan4 = msReaderParquet.getScanPoints(scanNumber + 3);
    QCOMPARE(scan4.first, eNoError);

    const QVector<ScanPoints*> scans = {
            scan1.second,
            scan2.second,
            scan3.second,
            scan4.second
    };

    MsScanInfo msScanInfo1;
    e = msReaderParquet.getMsScanInfo(scanNumber, &msScanInfo1);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo2;
    e = msReaderParquet.getMsScanInfo(scanNumber + 1, &msScanInfo2);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo3;
    e = msReaderParquet.getMsScanInfo(scanNumber + 2, &msScanInfo3);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo4;
    e = msReaderParquet.getMsScanInfo(scanNumber + 3, &msScanInfo4);
    QCOMPARE(e, eNoError);

    const QVector<MsScanInfo> msScanInfos = {
            msScanInfo1,
            msScanInfo2,
            msScanInfo3,
            msScanInfo4
    };


    DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
    e = deMultiplexScanerator.deMultiplexScans(scans, msScanInfos);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(DeMultiplexScaneratorTests)

#include "DeMultiplexScaneratorTests.moc"
