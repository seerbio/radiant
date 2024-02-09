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

//  TODO make this into a truncated file starting at scan 21796 to 21810
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

    QVector<QVector<float>> transitionMatrixVecs;
    e =  deMultiplexScanerator._buildTransitionMatrixTestAccess(scans, &transitionMatrixVecs);
    QCOMPARE(e, eNoError);

    const QVector<QVector<int>> expectedResults = {
            {10590, 120762, 656679, 0},
            {7947, 28843, 584602, 0},
            {0, 86665, 612637, 9280},
            {0, 1225137, 405833, 0}
    };

    QVector<QVector<int>> transitionMatrixVecsInts(4);
    for(int i = 0; i < transitionMatrixVecs.size(); i++) {
        for(int j = 0; j < transitionMatrixVecs.front().size(); j++) {
            transitionMatrixVecsInts[i].push_back(static_cast<int>(transitionMatrixVecs.at(i).at(j)));
        }
    }

    QCOMPARE(transitionMatrixVecsInts, expectedResults);

}

void DeMultiplexScaneratorTests::buildScanMaskMatrixTest() {

    ERR_INIT

    //  TODO make this into a truncated file starting at scan 21796 to 21810
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
            msScanInfo3,
            msScanInfo1,
            msScanInfo2,
            msScanInfo4
    };

    QVector<QVector<float>> scanMaskMatrixVecs;

    DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
    e = deMultiplexScanerator._buildScanMaskMatrixTestAccess(msScanInfos, &scanMaskMatrixVecs);
    QCOMPARE(e, eNoError);

    const QVector<QVector<float>> expectedResults = {
            {0, 0, 1, 1, 1, 1, 0},
            {1, 1, 1, 1, 0, 0, 0},
            {0, 1, 1, 1, 1, 0, 0},
            {0, 0, 0, 1, 1, 1, 1}
    };
    QCOMPARE(scanMaskMatrixVecs, expectedResults);


}

void DeMultiplexScaneratorTests::deMultiplexScansTest() {

    ERR_INIT

//  TODO make this into a truncated file starting at scan 21796 to 21810
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

    QVector<ScanPoints> demultiplexedScans;
    QVector<MzTargetKey> mzTargetKeys;
    double windwoSize;

    DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
    e = deMultiplexScanerator.deMultiplexScans(
            scans,
            msScanInfos,
            &demultiplexedScans,
            &mzTargetKeys,
            &windwoSize
            );
    QCOMPARE(e, eNoError);

    const QVector<MzTargetKey> mzTargetKeysExpected = {
            "501500",
            "502500",
            "503500",
            "504500",
            "505500",
            "506500",
            "507500"
    };

    QCOMPARE(demultiplexedScans.at(1).at(0), PointFF(145.052, 7947.07));
    QCOMPARE(mzTargetKeys, mzTargetKeysExpected);
    QCOMPARE(windwoSize, 1.0);

}


QTEST_MAIN(DeMultiplexScaneratorTests)

#include "DeMultiplexScaneratorTests.moc"
