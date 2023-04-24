#include "MsFrameScoretron.h"

#include <QElapsedTimer>
#include <QtTest>


class MsFrameScoretronTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronTests() = default;
    ~MsFrameScoretronTests() override = default;


private slots:

    void scoreCandidatesTest();

private:

    static PythiaParameters pythiaParameters() {


        PythiaParameters pythiaParameters;

        pythiaParameters.returnPSMTopN = 1;
        pythiaParameters.maxTandemPointCount = 2;
        pythiaParameters.ms2ExtractionWidthPPM = 12.0;
        pythiaParameters.precursorExtractionWindowThomsons = 1.0;
        pythiaParameters.chargeStateMin = 2;
        pythiaParameters.chargeStateMax = 3;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
        );

        return pythiaParameters;
    }



};


void MsFrameScoretronTests::scoreCandidatesTest() {

    ERR_INIT

    const QString mzMLFileURI
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib");


    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "474966";
    QPair<Err, QPair<UniqueMsInfoScanKey, QString>> result = MsFrameScoretron::scoreCandidates(
            pythiaParameters(),
            mzMLFileURI,
            fragLibPath,
            uniqueMsInfoScanKey,
            {474.966 - 5.5, 474.966 + 5.5}
            );


}


QTEST_MAIN(MsFrameScoretronTests)

#include "MsFrameScoretronTests.moc"
