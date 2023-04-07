#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "MsScansDenoiseTron.h"
#include "MsReaderMzML.h"


using namespace Error;


class MsScansDenoiseTronTests : public QObject
{
    Q_OBJECT

public:

    MsScansDenoiseTronTests() = default;
    ~MsScansDenoiseTronTests() override = default;

private Q_SLOTS:

    void denoiseMsScansTest();

private:

    //TODO use proper path procedures.
    const QString m_filepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

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

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
        );

        return pythiaParameters;
    }


};

void MsScansDenoiseTronTests::denoiseMsScansTest() {

    ERR_INIT

//    MsReaderMzML reader;
//    e = reader.openFile(m_filepath);
//    QCOMPARE(e, Error::eNoError);
//
//    QVector<TandemScanIon> tandemScanIons;
//    e = reader.tandemScanIons(&tandemScanIons);
//    QCOMPARE(e, Error::eNoError);
//
//    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTandemFrames;
//    e = MsReaderBase::sortDIATandemScansByMzTarget(
//            tandemScanIons,
//            &diaTandemFrames
//            );
//    QCOMPARE(e, Error::eNoError);
//
//    QCOMPARE(diaTandemFrames.size(), 62);
//
//    MsScansDenoiseTron denoiser;
//    e = denoiser.init(pythiaParameters());
//    QCOMPARE(e, Error::eNoError);
//
//    QMap<ScanNumber, ScanPoints> scansFrameDenoised;
//    e = denoiser.denoiseScansFrame(
//            diaTandemFrames.first(),
//            &scansFrameDenoised
//            );
//    QCOMPARE(e, Error::eNoError);
//
//    QCOMPARE(scansFrameDenoised.size(), 413);
}


QTEST_MAIN(MsScansDenoiseTronTests)

#include "MsScansDenoiseTronTests.moc"
