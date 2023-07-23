#include "IRTPredictron.h"

#include <QtTest/QtTest>


class IRTPredictronTests : public QObject
{
    Q_OBJECT

public:
    IRTPredictronTests() = default;
    ~IRTPredictronTests() = default;

private Q_SLOTS:
    void loadModelTest();
    void testPrediction();
    void testPredictionFail();

};

void IRTPredictronTests::loadModelTest()
{

    ERR_INIT

    const QString testModel = QDir(qApp->applicationDirPath()).filePath("iRT_Model.json");
    const QString testModelFail = QStringLiteral("KalliopeAndBellatrix.json");

    IRTPredictron predictotron;

    e = predictotron.init(testModelFail);
    QCOMPARE(e, eFileError);

    e = predictotron.init(testModel);
    QCOMPARE(e, eNoError);

}

void IRTPredictronTests::testPrediction() {

    ERR_INIT

    const QStringList peptideStringWithModsList = {
            "ASDLK",
            "ASDLGK",
    };

    const QString testModel = QDir(qApp->applicationDirPath()).filePath("iRT_Model.json");

    IRTPredictron predictotron;
    e = predictotron.init(testModel);
    QCOMPARE(e, eNoError);

    QVector<float> rawPredictionResults;
    e = predictotron.batchPredictIRT(peptideStringWithModsList, &rawPredictionResults);
    QCOMPARE(e, eNoError);

    qDebug() << rawPredictionResults;

}

void IRTPredictronTests::testPredictionFail() {

    ERR_INIT

    const QStringList peptideStringWithModsList = {
            "ASDLK",
            "EASDLGEDKASDLGEDKASDLGEDKASDLGEDKASDLGEDKASDLGEDK",
    };

    const QString testModel = QDir(qApp->applicationDirPath()).filePath("iRT_Model.json");

    IRTPredictron predictotron;
    e = predictotron.init(testModel);
    QCOMPARE(e, eNoError);

    QVector<float> rawPredictionResults;
    e = predictotron.batchPredictIRT(peptideStringWithModsList, &rawPredictionResults);
    QCOMPARE(e, eError);

}


QTEST_MAIN(IRTPredictronTests)
#include "IRTPredictronTests.moc"
