#include "TargetDecoyCandidatePairScoretron.h"

#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"

#include <QtTest/QtTest>

#include <iostream>
#include <set>

class TargetDecoyCandidatePairScoretronTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairScoretronTests() = default;
    ~TargetDecoyCandidatePairScoretronTests() override = default;

private Q_SLOTS:
    void calculateChunkSizingTest();
    void calculateSliceCountTest();
    void calculateSliceBoundsCoverageTest();
    void loadModelTest();


};

void TargetDecoyCandidatePairScoretronTests::calculateChunkSizingTest() {
    const int targetChunkSize
        = TargetDecoyCandidatePairScoretronUtils::calculateTargetChunkSize(100, 5);
    QCOMPARE(targetChunkSize, 5);

    const int minChunkSize
        = TargetDecoyCandidatePairScoretronUtils::calculateMinChunkSize(targetChunkSize);
    QCOMPARE(minChunkSize, 3);
}

void TargetDecoyCandidatePairScoretronTests::calculateSliceCountTest() {
    QCOMPARE(
        TargetDecoyCandidatePairScoretronUtils::calculateSliceCount(2, 5, 3),
        1
        );
    QCOMPARE(
        TargetDecoyCandidatePairScoretronUtils::calculateSliceCount(20, 5, 3),
        4
        );
    QCOMPARE(
        TargetDecoyCandidatePairScoretronUtils::calculateSliceCount(11, 3, 2),
        4
        );
}

void TargetDecoyCandidatePairScoretronTests::calculateSliceBoundsCoverageTest() {
    std::set<int> coveredIndexes;

    for (int sliceIndex = 0; sliceIndex < 4; ++sliceIndex) {
        const QPair<int, int> sliceBounds
            = TargetDecoyCandidatePairScoretronUtils::calculateSliceBounds(10, sliceIndex, 4);
        for (int itemIndex = sliceBounds.first;
             itemIndex < sliceBounds.first + sliceBounds.second;
             ++itemIndex) {
            const bool inserted = coveredIndexes.insert(itemIndex).second;
            QVERIFY(inserted);
        }
    }

    QCOMPARE(static_cast<int>(coveredIndexes.size()), 10);
    for (int itemIndex = 0; itemIndex < 10; ++itemIndex) {
        QVERIFY(coveredIndexes.find(itemIndex) != coveredIndexes.end());
    }

    const QPair<int, int> firstOddSlice
        = TargetDecoyCandidatePairScoretronUtils::calculateSliceBounds(5, 0, 2);
    const QPair<int, int> secondOddSlice
        = TargetDecoyCandidatePairScoretronUtils::calculateSliceBounds(5, 1, 2);
    QCOMPARE(firstOddSlice.first, 0);
    QCOMPARE(firstOddSlice.second, 2);
    QCOMPARE(secondOddSlice.first, 2);
    QCOMPARE(secondOddSlice.second, 3);
}

void TargetDecoyCandidatePairScoretronTests::loadModelTest() {

    QSKIP("TODO: use bundled test data");

    ERR_INIT

    QSKIP("reactivate when file pathing is figured out.");

    const QString fragLibUri
        = QStringLiteral("/home/anichols/Downloads/human_plasma_arath_entrapment.fasta.predicted.speclib.fragLib");

//    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
//    e = targetDecoyCandidatePairManager.init(
//            PythiaParameterReader::genericPythiaParametersForTests(),
//            fragLibUri
//            );
//    QCOMPARE(e, eNoError);
//
//    const QString msDataFilePath
//        = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.prq");
////        = QStringLiteral("/home/anichols/Downloads/EXP23109_2023astral006cX26_A.raw.mzML");
//
//    MsReaderPointerAcc msReaderPointerAcc;
//    e = msReaderPointerAcc.openFile(msDataFilePath);
//    QCOMPARE(e, eNoError);
//
//    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrame;
//    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(
//            &diaTargetFrame
//    );
//    QCOMPARE(e, eNoError);
//
//    const int msLevel = 1;
//    QMap<ScanNumber, ScanPoints*> scanNumberVsScanTimeMS1;
//    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanTimeMS1);
//    QCOMPARE(e, eNoError);
//
//    TargetDecoyCandidatePairScoretron targetDecoyCandidatePairScoretron;
//    e = targetDecoyCandidatePairScoretron.init(
//            PythiaParameterReader::genericPythiaParametersForTests(),
//            scanNumberVsScanTimeMS1,
//            &msReaderPointerAcc,
//            &diaTargetFrame,
//            &targetDecoyCandidatePairManager
//            );
//    QCOMPARE(e, eNoError);
//
//    MsCalibratomatic msCalibratomatic;
//
//    int topNMS2Ions = 6;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
//    e = targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
//            topNMS2Ions,
//            msCalibratomatic,
//            &scoredTargetDecoyPointers
//            );
//    QCOMPARE(e, eNoError);
//    QCOMPARE(scoredTargetDecoyPointers.size(), 228217);

}


QTEST_MAIN(TargetDecoyCandidatePairScoretronTests)
#include "TargetDecoyCandidatePairScoretronTests.moc"
