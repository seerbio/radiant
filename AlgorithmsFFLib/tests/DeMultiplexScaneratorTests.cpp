#include <QtTest>
#include <QCoreApplication>

#include "DeMultiplexScanerator.h"
#include "MsReaderParquet.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

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

//    static void deMultiplexEntireFileResearch();

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
//
//namespace {
//
//    struct DemultiParallelInput {
//        QVector<ScanPoints*> scans;
//        QVector<MsScanInfo> msScanInfos;
//    };
//
//    struct DemultiParallelOutput {
//        QVector<ScanPoints> demultiplexedScans;
//        QVector<MzTargetKey> mzTargetKeys;
//        double windwoSize;
//    };
//
//    QPair<Err, DemultiParallelOutput> demultiplexLogic(const DemultiParallelInput &demultiParallelInput) {
//
//        ERR_INIT
//
//        DemultiParallelOutput demultiParallelOutput;
//        DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
//        e = deMultiplexScanerator.deMultiplexScans(
//                demultiParallelInput.scans,
//                demultiParallelInput.msScanInfos,
//                &demultiParallelOutput.demultiplexedScans,
//                &demultiParallelOutput.mzTargetKeys,
//                &demultiParallelOutput.windwoSize
//        );
//
//        return {e, demultiParallelOutput};
//    }
//
//}//namespace
//void DeMultiplexScaneratorTests::deMultiplexEntireFileResearch() {
//
//    ERR_INIT
//
//    //  TODO make this into a truncated file starting at scan 21796 to 21810
//    const QString prqFFFile = "/home/anichols/Downloads/TESTING_20240130RC6_30min_iso-4-500-510-slide-7test-15k.mzML.prqFF";
//
//    MsReaderParquet msReaderParquet;
//    e = msReaderParquet.openFile(prqFFFile);
//    QCOMPARE(e, eNoError);
//
//    const int msLevel = 2;
//    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderParquet.getMsScanInfos(msLevel);
//
//    const QHash<MzTargetKey, bool> lastTargetKeys = {
//            {"506000", true}
//    };
//
//    QVector<QVector<MsScanInfo>> groupedMsScanInfosAll;
//
//    QVector<MsScanInfo> groupedMsScanInfosCurrent;
//    for(const MsScanInfo &msi : msScanInfos) {
//        if (!lastTargetKeys.contains(msi.targetKey())) {
//            groupedMsScanInfosCurrent.push_back(msi);
//            continue;
//        }
//
//        groupedMsScanInfosCurrent.push_back(msi);
//        groupedMsScanInfosAll.push_back(groupedMsScanInfosCurrent);
//        groupedMsScanInfosCurrent.clear();
//    }
//
//    QMap<ScanNumber, ScanPoints*> scanPoints;
//    e = msReaderParquet.getScanPoints(2, &scanPoints);
//    QCOMPARE(e, eNoError);
//
//    QVector<QVector<ScanPoints*>> groupedScanPointsAll;
//    QVector<ScanPoints*> groupedScanPointsCurrent;
//    for (const QVector<MsScanInfo> &msis : groupedMsScanInfosAll) {
//
//        for (const MsScanInfo &m : msis) {
//            ScanPoints *sp = scanPoints.value(m.scanNumber);
//            groupedScanPointsCurrent.push_back(sp);
//
//        }
//
//        groupedScanPointsAll.push_back(groupedScanPointsCurrent);
//        groupedScanPointsCurrent.clear();
//    }
//
//    e = ErrorUtils::isEqual(groupedMsScanInfosAll.size(), groupedScanPointsAll.size());
//    QCOMPARE(e, eNoError);
//
//    QVector<DemultiParallelInput> inputs;
//    for (int i = 0; i < groupedScanPointsAll.size(); i++) {
//        DemultiParallelInput demultiParallelInput;
//        demultiParallelInput.msScanInfos = groupedMsScanInfosAll.at(i);
//        demultiParallelInput.scans = groupedScanPointsAll.at(i);
//
//        inputs.push_back(demultiParallelInput);
//    }
//
//    QFuture<QPair<Err, DemultiParallelOutput>> results = QtConcurrent::mapped(
//            inputs,
//            demultiplexLogic
//            );
//    results.waitForFinished();
//
//    QMap<ScanNumber, ScanPoints> newScanPoints;
//    QMap<ScanNumber, MsScanInfo> newMsScanInfos;
//
//    int newScanNumber = 0;
//    int outputCounter = 0;
//    for (const QPair<Err, DemultiParallelOutput> &res : results) {
//        e = res.first;
//        QCOMPARE(e, eNoError);
//
//        const DemultiParallelInput &demultiParallelInput = inputs.at(outputCounter++);
//        const DemultiParallelOutput &demultiParallelOutput = res.second;
//
//        e = ErrorUtils::isEqual(
//                demultiParallelOutput.mzTargetKeys.size(),
//                demultiParallelOutput.demultiplexedScans.size()
//                );
//        QCOMPARE(e, eNoError);
//
//        for (int i = 0; i < demultiParallelOutput.mzTargetKeys.size(); i++) {
//
//            MsScanInfo newMsScanInfo;
//            newMsScanInfo.isoWindowLower = demultiParallelOutput.windwoSize / 2.0;
//            newMsScanInfo.isoWindowUpper = demultiParallelOutput.windwoSize / 2.0;
//            newMsScanInfo.scanNumber = newScanNumber++;
//            newMsScanInfo.msLevel = 2;
//            newMsScanInfo.ionMobilityIndex = -1;
//            newMsScanInfo.ionMobilityDriftTime = -1.0;
//            newMsScanInfo.scanTime = demultiParallelInput.msScanInfos.front().scanTime;
//            newMsScanInfo.collisionEnergy = demultiParallelInput.msScanInfos.front().collisionEnergy;
//            e = ErrorUtils::toFloat(demultiParallelOutput.mzTargetKeys.at(i) , &newMsScanInfo.precursorTargetMz);
//            QCOMPARE(e, eNoError);
//
//            newMsScanInfos.insert(newMsScanInfo.scanNumber, newMsScanInfo);
//            newScanPoints.insert(newMsScanInfo.scanNumber, demultiParallelOutput.demultiplexedScans.at(i));
//        }
//
//    }
//
//    msReaderParquet.setMsScanInfo(newMsScanInfos);
//    e = msReaderParquet.setScanPoints(newScanPoints);
//    QCOMPARE(e, eNoError);
//
//}


QTEST_MAIN(DeMultiplexScaneratorTests)

#include "DeMultiplexScaneratorTests.moc"
