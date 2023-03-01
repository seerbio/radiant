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

private:

    static PythiaParameters pythiaParameters() {

        const QString paramsFile = "/home/anichols/Repositories/PythiaDIACpp/FileReadersLib/tests/TestFiles/WorkFlowTestsParams.pythia";

        PythiaParameterReader reader;
        PythiaParameters pythiaParameters;
        reader.readFile(paramsFile);
        reader.loadPythiaParameters(&pythiaParameters);

        return pythiaParameters;
    }

};

void PythiaDIAWorkflowTests::execTest() {

    ERR_INIT

    const QString mzMLFileURI
        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.fragLib");

    const QString pepLibPath
            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.pepLib");

    PythiaDIAWorkflow pythiaDiaWorkflow;
    e = pythiaDiaWorkflow.init(
            pythiaParameters(),
            fragLibPath,
            pepLibPath
            );
    QCOMPARE(e, eNoError);

    e = pythiaDiaWorkflow.processFile(mzMLFileURI);
    QCOMPARE(e, eNoError);


}


QTEST_MAIN(PythiaDIAWorkflowTests)

#include "PythiaDIAWorkflowTests.moc"
