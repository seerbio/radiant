#include "Ms1FeatureFinder.h"

#include "FeatureFinderHill.h"
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

};

void Ms1FeatureFinderTests::execTest() {

    ERR_INIT

    const QString prqfilepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.reCal");
//            = "/home/anichols/Downloads/EXP20011_2023ms0452X12_B.raw.mzML.prq";

    Ms1FeatureFinder featureFinder;
    e = featureFinder.init(PythiaParameterReader::genericPythiaParametersForTests());
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHill> featureFinderHills;
    e = featureFinder.exec(prqfilepath, &featureFinderHills);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(Ms1FeatureFinderTests)

#include "Ms1FeatureFinderTests.moc"
