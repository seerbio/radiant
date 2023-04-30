#include "Ms1FeatureFinder.h"

#include "MsFrame.h"
#include "MsReaderParquet.h"

#include <QElapsedTimer>
#include <QVector>
#include <QtTest>


class Ms1FeatureFinderTests : public QObject
{

    Q_OBJECT
    
public:
    Ms1FeatureFinderTests() = default;
    ~Ms1FeatureFinderTests() override = default;


private slots:

    void execTest();

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


void Ms1FeatureFinderTests::execTest() {

    ERR_INIT

    const QString prqfilepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.reCal");


    Ms1FeatureFinder featureFinder;
    e = featureFinder.init(pythiaParameters());
    QCOMPARE(e, eNoError);

    e = featureFinder.exec(prqfilepath);
    QCOMPARE(e, eNoError);






}


QTEST_MAIN(Ms1FeatureFinderTests)

#include "Ms1FeatureFinderTests.moc"
