#include "MsFraggerTronWorkFlow.h"

#include "Error.h"
#include "MsFraggerTronResultsReader.h"
#include "MsReaderBase.h"
#include "PeptidesLibraryTron.h"
#include "PythiaParameterReader.h"

#include <QtTest>

class MsFraggerTronWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    MsFraggerTronWorkFlowTests() = default;
    ~MsFraggerTronWorkFlowTests() override = default;

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

        return pythiaParameters;
    }

};

void MsFraggerTronWorkFlowTests::execTest() {

    ERR_INIT

    const QString mzMLFileURI
        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/RawData/human_plasma_entrapment_super_trunc.fasta.fragLib");

    const QString pepLibPath
            = QStringLiteral("/home/anichols/Desktop/RawData/human_plasma_entrapment_super_trunc.fasta.pepLib");

    MsFraggerTronWorkFlow msFraggerTronWorkFlow;
    e = msFraggerTronWorkFlow.init(
            pythiaParameters(),
            fragLibPath,
            pepLibPath
            );
    QCOMPARE(e, eNoError);

    QString firstPassPSMsFilePath;
    QVector<TandemScanIon> tandemScanIons;
    e = msFraggerTronWorkFlow.processFile(
            mzMLFileURI,
            &firstPassPSMsFilePath,
            &tandemScanIons
            );
    QCOMPARE(e, eNoError);

    QVector<RowToWrite> rowsToWrite;
    e = MsFraggerTronResultsReader::readCsv(
            firstPassPSMsFilePath,
            &rowsToWrite
            );
    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsFraggerTronWorkFlowTests)

#include "MsFraggerTronWorkFlowTests.moc"
