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

        const QString &paramsFile
                = QDir(qApp->applicationDirPath()).filePath("WorkFlowTestsParams.pythia");

        PythiaParameterReader reader;
        PythiaParameters pythiaParameters;
        reader.readFile(paramsFile);
        reader.loadPythiaParameters(&pythiaParameters);

        pythiaParameters.topNMs2Ions = 12;
//        pythiaParameters.ms2ExtractionWidthPPM = 20;

        return pythiaParameters;
    }

};

void PythiaDIAWorkflowTests::execTest() {

    ERR_INIT

    const QString mzMLFileURI
        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.fragLib");
            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";

    const QString pepLibPath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/human_plasma_entrapment_super_trunc.fasta.pepLib");
            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.pepLib";

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
