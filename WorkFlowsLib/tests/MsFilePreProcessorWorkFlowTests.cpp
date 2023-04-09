#include "Error.h"
#include "MsFilePreProcessorWorkFlow.h"
#include "MsReaderPointerFactory.h"
#include "PythiaParameterReader.h"

#include <QtTest>


using namespace Error;


class MsFilePreProcessorWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    MsFilePreProcessorWorkFlowTests() = default;
    ~MsFilePreProcessorWorkFlowTests() override = default;

private slots:

    void preprocessTandemScansTest();

};

void MsFilePreProcessorWorkFlowTests::preprocessTandemScansTest() {

    ERR_INIT

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
    params.featureFinderTolerancePPM = 12;
    params.skipScanCount = 2;
    params.minScanCount = 3;
    params.useMeanMz = true;

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

    const QString &prqFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq";
    QPair<Err, MsReaderPointer> msReader = MsReaderPointerFactory::createInstance(prqFilePath);
    QCOMPARE(msReader.first, eNoError);

}


QTEST_MAIN(MsFilePreProcessorWorkFlowTests)

#include "MsFilePreProcessorWorkFlowTests.moc"
