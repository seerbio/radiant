#include "FDRCLassifierNeuralNet.h"
#include "CandidateScores.h"

#include <QtTest/QtTest>

#include <iostream>

class FDRCLassifierNeuralNetTests : public QObject
{
    Q_OBJECT

public:
    FDRCLassifierNeuralNetTests() = default;
    ~FDRCLassifierNeuralNetTests() override = default;

private Q_SLOTS:
    void initTest();
    void execTest();
};

void FDRCLassifierNeuralNetTests::initTest() {

    ERR_INIT

    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = fdrClassifierNeuralNet.init(
            5,
            1,
            5,
            0.1,
            0.01,
            {1,20}
            );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            0,
            5,
            0.1,
            0.01,
            {1,20}
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            0,
            0.1,
            0.01,
            {1,20}
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.0,
            0.01,
            {1,20}
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.1,
            0.0,
            {1,20}
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.1,
            0.0000001,
            {1,20}
    );
    QCOMPARE(e, eNoError);

}

void FDRCLassifierNeuralNetTests::execTest() {
    QSKIP("TODO: enable with internal test data");

    ERR_INIT
//TODO figure out proper filepathing for tests.
//    const QString &scoredCandsFullFragsFilePath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/scoredCandidatesAllFullFragIons.parquet");
//
//    const QString &keyVsScoredCandidateCulledFilePath
//            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/keyVsScoredCandidateCulled.parquet");
//
//    QVector<CandidateScores> scoredCandidatesAllFullFragIons;
//    e = ParquetReader::read(scoredCandsFullFragsFilePath, &scoredCandidatesAllFullFragIons);
//    QCOMPARE(e, eNoError);
//
//    QVector<CandidateScores> scoredCandidateCulledVec;
//    e = ParquetReader::read(keyVsScoredCandidateCulledFilePath, &scoredCandidateCulledVec);
//    QCOMPARE(e, eNoError);
//
//    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
//    e = fdrClassifierNeuralNet.init(
//            12,
//            10,
//            6,
//            0.1,
//            0.001,
//            {1,20}
//    );
//    QCOMPARE(e, eNoError);
//
//    QVector<CandidateScores> scoredCandidatesClassifier;
//    e = fdrClassifierNeuralNet.exec(
//            scoredCandidateCulledVec,
//            &scoredCandidatesClassifier
//            );
//    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FDRCLassifierNeuralNetTests)
#include "FDRCLassifierNeuralNetTests.moc"
