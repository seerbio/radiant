#include "TargetDecoyCandidatePairScoretron.h"

#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"

#include <QtTest/QtTest>

#include <iostream>

class TargetDecoyCandidatePairScoretronTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairScoretronTests() = default;
    ~TargetDecoyCandidatePairScoretronTests() override = default;

private Q_SLOTS:
    void loadModelTest();


};

void TargetDecoyCandidatePairScoretronTests::loadModelTest() {

    ERR_INIT

    const QString fragLibUri
        = QStringLiteral("/home/anichols/Downloads/human_plasma_arath_entrapment.fasta.predicted.speclib.fragLib");

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            fragLibUri
            );
    QCOMPARE(e, eNoError);

    const QString msDataFilePath
        = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrame;
    e = msReaderPointerAcc.ptr->collateTandemPrecursorTargetsDIA(
            &diaTargetFrame
    );
    QCOMPARE(e, eNoError);

    TargetDecoyCandidatePairScoretron targetDecoyCandidatePairScoretron;
    e = targetDecoyCandidatePairScoretron.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            &msReaderPointerAcc,
            &diaTargetFrame,
            &targetDecoyCandidatePairManager
            );
    QCOMPARE(e, eNoError);

    MsCalibratomatic msCalibratomatic;

    int topNMS2Ions = 6;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
    e = targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMS2Ions,
            0.2,
            msCalibratomatic,
            &scoredTargetDecoyPointers
            );
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(TargetDecoyCandidatePairScoretronTests)
#include "TargetDecoyCandidatePairScoretronTests.moc"
