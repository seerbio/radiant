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
//            = "/home/anichols/Desktop/RawData/2022-05-05-decoys-Uniprot_human_plus_Arabidopsis.fasta";
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

    Modification carboxyAmidoMethyl(
            'C',
            "CAM",
            ModificationType::FIXED,
            "C2H3NO"
            );

    params.modifications = {carboxyAmidoMethyl};
    e = PythiaParameterReader::applyFixedModificationsToAminoAcids(
            params,
            &params.aminoAcids
            );
    QCOMPARE(e, eNoError);

    FastaFileToPeptidesListWorkFlow wf;
    e = wf.init(params);
    QCOMPARE(e, eNoError);

    QString outputFilePath;
    e = wf.exec(
            fastaFilePath,
            &outputFilePath
            );
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FastaFileToPeptidesListWorkFlowTests)

#include "FastaFileToPeptidesListWorkFlowTests.moc"
