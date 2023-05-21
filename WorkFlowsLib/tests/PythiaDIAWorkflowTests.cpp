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


};

void PythiaDIAWorkflowTests::execTest() {

    ERR_INIT

    const QString mzMLFileURI
//        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");
        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.reCal");

    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.fragLib");
//            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";
            = "/home/anichols/Desktop/2022_02_22_Homo_sapiens_UP000005640.fragLib";

    PythiaDIAWorkflow pythiaDiaWorkflow;
    e = pythiaDiaWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            fragLibPath
            );
    QCOMPARE(e, eNoError);

    e = pythiaDiaWorkflow.processFile(mzMLFileURI);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(PythiaDIAWorkflowTests)

#include "PythiaDIAWorkflowTests.moc"
