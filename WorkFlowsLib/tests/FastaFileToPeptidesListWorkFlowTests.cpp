#include "Error.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "PythiaParameterReader.h"

#include <QtTest>


using namespace Error;


class FastaFileToPeptidesListWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    FastaFileToPeptidesListWorkFlowTests() = default;
    ~FastaFileToPeptidesListWorkFlowTests() override = default;

private slots:

    void execTest();

};

void FastaFileToPeptidesListWorkFlowTests::execTest() {

    ERR_INIT

    const QString &fastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    const QString &fastaLarge
        = "/home/anichols/Desktop/RawData/UP000000589_10090.target.decoy.con.pepcal.20211006.fasta";


    PythiaParameters params;
    params.cTermCleavePoints = QStringList({"K", "R"});
    params.addDecoys = false;
    params.peptideLengthMin = 7;
    params.peptideLengthMax = 40;
    params.chargeStateMin = 2;
    params.chargeStateMax = 3;
    params.maxTandemPointCount = 500;
    params.returnPSMTopN = 1;
    params.ms2ExtractionWidthPPM = 12.0;
    params.precursorExtractionWindowThomsons = 1.0;

    FastaFileToPeptidesListWorkFlow wf;
    e = wf.init(params);
    QCOMPARE(e, eNoError);

    QString outputFilePath;
    e = wf.exec(
            fastaLarge,
            &outputFilePath
            );
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FastaFileToPeptidesListWorkFlowTests)

#include "FastaFileToPeptidesListWorkFlowTests.moc"
