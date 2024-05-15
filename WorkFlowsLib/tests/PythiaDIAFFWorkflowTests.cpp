#include "PythiaDIAFFWorkflow.h"

#include "Error.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

#include <QtTest>

class PythiaDIAFFWorkflowTests : public QObject
{
    Q_OBJECT

public:

    PythiaDIAFFWorkflowTests() = default;
    ~PythiaDIAFFWorkflowTests() override = default;

private slots:

    void initTest();

    void buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest();


};

void PythiaDIAFFWorkflowTests::initTest() {

    ERR_INIT

    const QString &testFragLibFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    const QString &testFastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    PythiaDIAFFWorkflow pythiaDiaffWorkflow;
    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            testFastaFilePath
            );
    QCOMPARE(e, eNoError);

    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            "kalliope.fragLibDF",
            testFastaFilePath
    );
    QCOMPARE(e, eFileError);

    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            "bellatrix.fasta"
    );
    QCOMPARE(e, eFileError);

}

void PythiaDIAFFWorkflowTests::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest() {

    ERR_INIT

    const QString &testFragLibFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    const QString &testFastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    const QString &testMsFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(testMsFilePath);
    QCOMPARE(e, eNoError);

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();

    PythiaDIAFFWorkflow pythiaDiaffWorkflow;
    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            testFastaFilePath
            );
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    // e = pythiaDiaffWorkflow.buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
    //         uniqueMsScanInfos,
    //         -1.0,
    //         &mzTargetKeyVsTargetDecoyCandidatePointers
    //         );
    QCOMPARE(e, eNoError);
    QCOMPARE(mzTargetKeyVsTargetDecoyCandidatePointers.value("994850").size(), 65);

    QVector<TargetDecoyCandidatePair*> pairs = mzTargetKeyVsTargetDecoyCandidatePointers.value("994850");

    std::sort(
            pairs.begin(),
            pairs.end(),
            [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){return l->mz() < r->mz();}
            );

    PeptideStringWithMods peptideStringWithMods = pairs.at(0)->peptideStringWithMods();
    QCOMPARE(peptideStringWithMods, "ASEAEDASLLSFMQGYMK");

}


QTEST_MAIN(PythiaDIAFFWorkflowTests)

#include "PythiaDIAFFWorkflowTests.moc"
