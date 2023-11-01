#include "IRTPredictron.h"
#include "MathUtils.h"
#include "NearestNeighborsSearch.h"

#include <QtTest/QtTest>

#include <iostream>

class IRTPredictronTests : public QObject
{
    Q_OBJECT

public:
    IRTPredictronTests() = default;
    ~IRTPredictronTests() override = default;

private Q_SLOTS:
    void loadModelTest();
    void testPrediction();
    void testPredictionFail();
    void buildNearestNeighborsIRTDataTest();
    void iRTToScanTimeTest();

};

void IRTPredictronTests::loadModelTest() {

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

    QVERIFY(MathUtils::tSame(rawPredictionResults.front(), static_cast<float>(14.3491)));
    QVERIFY(MathUtils::tSame(rawPredictionResults.back(), static_cast<float>(11.9481)));

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

void IRTPredictronTests::buildNearestNeighborsIRTDataTest() {

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    //TODO ad this to testing files and properly make the path.
    const QString iRTReCalFilePath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.iRT");

    QVector<QPair<double, Coors>> nnInputData;
//    e = IRTPredictron::buildNearestNeighborsIRTData(iRTReCalFilePath, &nnInputData);
//    QCOMPARE(e, eNoError);
//    QCOMPARE(2957, nnInputData.size());


}

void IRTPredictronTests::iRTToScanTimeTest() {

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    //TODO ad this to testing files and properly make the path.
    const QString iRTReCalFilePath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.iRT");

    QVector<QPair<double, Coors>> nnInputData;
//    e = IRTPredictron::buildNearestNeighborsIRTData(iRTReCalFilePath, &nnInputData);
//    QCOMPARE(e, eNoError);
//
//    const int midPoint = static_cast<int>(nnInputData.size() / 2.0);
//    QVector<QPair<double, Coors>> nnInputDataTrain = nnInputData.mid(0, midPoint);
//    QVector<QPair<double, Coors>> nnInputDataTest = nnInputData.mid(midPoint, midPoint);
//
//    qDebug() << nnInputDataTest.size() << nnInputDataTrain.size();
//
//    NearestNeighborsSearch nearestNeighborsSearch;
//    e = nearestNeighborsSearch.init(nnInputDataTrain);
//    QCOMPARE(e, eNoError);
//
//    const int kNearestPoints = 10;
//
//    for (const QPair<double, Coors> &row : nnInputDataTest) {
//
//        const double scanTime = row.first;
//        const Coors &coor = row.second;
//
//        QVector<NNSearchResult> nnSearchResult;
//        nearestNeighborsSearch.kNearestNeighborsSearch(
//                {coor},
//                kNearestPoints,
//                &nnSearchResult
//        );
//
////        std::cout << "(" << scanTime << "," << nnSearchResult.front().meanValues << "),"<< std::endl;
//
//    }

}


QTEST_MAIN(IRTPredictronTests)
#include "IRTPredictronTests.moc"
