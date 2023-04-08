#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "MsScansDenoiseTron.h"
#include "MsReaderParquet.h"


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
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

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

    MsReaderParquet reader;
    e = reader.openFile(m_filepath);
    QCOMPARE(e, Error::eNoError);

    QMap<ScanNumber, ScanPoints> ms1Scans;

    const int msLevel = 1;
    e = reader.getScanPoints(
            msLevel,
            &ms1Scans
    );
    QCOMPARE(e, eNoError);

    MsScansDenoiseTron denoiser;
    e = denoiser.init(pythiaParameters());
    QCOMPARE(e, Error::eNoError);

    QMap<ScanNumber, ScanPoints> scansFrameDenoised;
    e = denoiser.denoiseScansFrame(
            ms1Scans,
            &scansFrameDenoised
    );
    QCOMPARE(e, Error::eNoError);

    QCOMPARE(scansFrameDenoised.size(), 413);

    const QVector<int> expectedDenoisedScanSizePerScanNumber = {
        34, 209, 213, 208, 225, 185, 178, 188, 180, 196, 181, 189, 180, 182, 185, 190, 193, 203, 200,
        197, 201, 201, 204, 204, 193, 195, 198, 189, 191, 195, 185, 193, 174, 176, 200, 185, 186, 191,
        201, 180, 173, 176, 180, 173, 180, 181, 190, 175, 186, 183, 184, 183, 183, 191, 182, 172, 191, 180,
        186, 188, 175, 182, 206, 185, 203, 184, 194, 190, 188, 189, 189, 194, 189, 198, 193, 191, 226, 222,
        208, 225, 227, 226, 248, 226, 211, 198, 136, 179, 241, 243, 128, 143, 110, 79, 79, 125, 90, 120, 172,
        166, 160, 167, 201, 135, 152, 169, 239, 270, 332, 239, 237, 182, 229, 208, 188, 168, 196, 212, 243,
        270, 242, 269, 214, 224, 239, 290, 248, 344, 333, 316, 264, 262, 236, 206, 207, 208, 183, 230, 249,
        231, 251, 265, 280, 318, 266, 303, 259, 265, 240, 164, 217, 217, 222, 254, 225, 223, 110, 155, 259,
        271, 297, 294, 285, 311, 337, 271, 257, 264, 221, 294, 331, 374, 405, 426, 438, 428, 362, 286, 227,
        179, 231, 311, 323, 362, 380, 401, 365, 330, 302, 301, 320, 366, 360, 361, 440, 397, 343, 304, 376,
        414, 373, 336, 293, 277, 355, 365, 384, 307, 349, 339, 307, 325, 320, 336, 404, 500, 477, 500, 337,
        313, 242, 261, 276, 270, 365, 394, 524, 478, 541, 563, 527, 447, 305, 281, 251, 236, 234, 242, 276,
        341, 440, 483, 481, 420, 380, 431, 427, 422, 243, 264, 271, 335, 434, 487, 538, 490, 440, 447, 432,
        391, 417, 327, 391, 400, 446, 471, 551, 522, 401, 409, 325, 440, 351, 355, 357, 372, 201, 242, 314, 373,
        487, 581, 564, 477, 559, 532, 531, 539, 519, 348, 400, 365, 441, 509, 549, 547, 376, 377, 481, 498, 477,
        322, 293, 353, 309, 284, 362, 346, 328, 343, 377, 473, 535, 458, 493, 528, 404, 454, 411, 319, 202, 226,
        284, 307, 456, 566, 595, 580, 528, 297, 374, 444, 480, 424, 490, 511, 428, 394, 134, 139, 214, 338, 448,
        448, 464, 491, 391, 371, 309, 179, 253, 388, 211, 334, 544, 524, 603, 601, 535, 267, 380, 552, 390, 374,
        391, 497, 517, 602, 526, 513, 445, 463, 303, 344, 505, 517, 544, 519, 659, 662, 712, 719, 825, 859, 606,
        811, 804, 770, 878, 846, 925, 977, 1005, 1092, 907, 883, 930, 900, 842, 544, 281, 521, 388, 368, 386, 368,
        250, 389, 929, 1012, 1051, 948, 714
    };

    int counter = 0;
    for (const ScanNumber scanNumber : scansFrameDenoised.keys()) {
        QCOMPARE(scansFrameDenoised.value(scanNumber).size(), expectedDenoisedScanSizePerScanNumber.at(counter++));
    }
}


QTEST_MAIN(MsScansDenoiseTronTests)

#include "MsScansDenoiseTronTests.moc"
