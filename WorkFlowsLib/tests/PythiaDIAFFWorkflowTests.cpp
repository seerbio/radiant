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

}


QTEST_MAIN(PythiaDIAFFWorkflowTests)

#include "PythiaDIAFFWorkflowTests.moc"
