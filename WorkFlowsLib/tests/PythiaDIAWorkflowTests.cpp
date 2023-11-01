#include "PythiaDIAWorkflow.h"

#include "Error.h"
#include "PythiaParameterReader.h"

#include <QtTest>

class PythiaDIAWorkflowTests : public QObject
{
    Q_OBJECT

public:

    PythiaDIAWorkflowTests() = default;
    ~PythiaDIAWorkflowTests() override = default;

private slots:

    void execTest();
    void execIRTTest();


};

void PythiaDIAWorkflowTests::execTest() {

    QSKIP("TODO: enable with internal test data");

    ERR_INIT

    const QString mzMLFileURI
//        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");
        = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");

    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.fragLib");
//            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";
            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");

    const QString fragLibBackgroundPath
            = QStringLiteral("/home/anichols/Desktop/Libraries/uniparc_upid_UP000027126_2023_07_02.fasta.csv.fragLib");

//    PythiaDIAWorkflow pythiaDiaWorkflow;
//    e = pythiaDiaWorkflow.init(
//            PythiaParameterReader::genericPythiaParametersForTests(),
//            fragLibPath,
//            fragLibBackgroundPath
//            );
//    QCOMPARE(e, eNoError);
//
//    e = pythiaDiaWorkflow.processFile(mzMLFileURI);
//    QCOMPARE(e, eNoError);

}

void PythiaDIAWorkflowTests::execIRTTest() {

    QSKIP("TODO: enable with internal test data");

    ERR_INIT

    QSKIP("activate when proper pathing is used");

    const QString mzMLFileURI
//        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");

    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.fragLib");
//            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";
            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");

    const QString fragLibBackgroundPath
            = QStringLiteral("/home/anichols/Desktop/Libraries/uniparc_upid_UP000027126_2023_07_02.fasta.csv.fragLib");

    const QString iRTReCalFilePath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.iRT");

//    PythiaDIAWorkflow pythiaDiaWorkflow;
//    e = pythiaDiaWorkflow.init(
//            PythiaParameterReader::genericPythiaParametersForTests(),
//            fragLibPath,
//            iRTReCalFilePath
//    );
//    QCOMPARE(e, eNoError);
//
//    e = pythiaDiaWorkflow.processFile(mzMLFileURI);
//    QCOMPARE(e, eNoError);

}


QTEST_MAIN(PythiaDIAWorkflowTests)

#include "PythiaDIAWorkflowTests.moc"
