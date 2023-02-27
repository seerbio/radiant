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

        Modification modOx('M', "Ox", ModificationType::DYNAMIC, "O");
        Modification modDeam('N', "Deam", ModificationType::DYNAMIC, "N");
        Modification modCAM('C', "CAM", ModificationType::FIXED, "C2H3NO");
        Modification modAce("N-term-protein", "Ace", ModificationType::DYNAMIC, "C2H5O");

        PythiaParameters pythiaParameters;
        pythiaParameters.modifications.append({modCAM, modOx, modDeam, modAce});
        pythiaParameters.cTermCleavePoints.append({"K", "R"});
        pythiaParameters.allowedMissedCleavages = 1;
        pythiaParameters.addDecoys = true;
        pythiaParameters.maxModificationsPeptide = 2;
        pythiaParameters.ms2ExtractionWidthPPM = 12.0;
        pythiaParameters.precursorExtractionWindowThomsons = 1.0;
        pythiaParameters.chargeStateMin = 2;
        pythiaParameters.chargeStateMax = 3;
        pythiaParameters.returnPSMTopN = 2;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
                );

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

    e = msFraggerTronWorkFlow.processFile(mzMLFileURI);
    QCOMPARE(e, eNoError);

    const QString firstPassPSMsFilePath
        = mzMLFileURI + S_GLOBAL_SETTINGS.DOT_PSM + S_GLOBAL_SETTINGS.DOT_CSV;

    QVector<RowToWrite> rowsToWrite;
    MsFraggerTronResultsReader::readCsv(
            firstPassPSMsFilePath,
            &rowsToWrite
            );

}


QTEST_MAIN(MsFraggerTronWorkFlowTests)

#include "MsFraggerTronWorkFlowTests.moc"
