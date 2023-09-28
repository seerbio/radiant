#include "FDRCLassifierNeuralNet.h"

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
            0.01
            );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            0,
            5,
            0.1,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            0,
            0.1,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.0,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.1,
            0.0
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            6,
            1,
            1,
            0.1,
            0.0000001
    );
    QCOMPARE(e, eNoError);

}

void FDRCLassifierNeuralNetTests::execTest() {

    ERR_INIT

    const QString &scoredCandsFullFragsFilePath
            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/scoredCandidatesAllFullFragIons.parquet");

    const QString &keyVsScoredCandidateCulledFilePath
            = QStringLiteral("/home/anichols/Repositories/Builds/PythiaDIACpp/bin/keyVsScoredCandidateCulled.parquet");

    QVector<ScoredCandidate> scoredCandidatesAllFullFragIons;
    e = ParquetReader::read(scoredCandsFullFragsFilePath, &scoredCandidatesAllFullFragIons);
    QCOMPARE(e, eNoError);

    QVector<ScoredCandidate> keyVsScoredCandidateCulledVec;
    e = ParquetReader::read(keyVsScoredCandidateCulledFilePath, &keyVsScoredCandidateCulledVec);
    QCOMPARE(e, eNoError);

    QMap<QString, ScoredCandidate> keyVsScoredCandidateCulled;
    for (const ScoredCandidate &sc : keyVsScoredCandidateCulledVec) {
        const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                sc.peptideStringWithMods,
                sc.targetKey,
                sc.charge
        );
        keyVsScoredCandidateCulled.insert(key, sc);
    }

    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = fdrClassifierNeuralNet.init(
            12,
            10,
            6,
            0.1,
            0.001
    );
    QCOMPARE(e, eNoError);

    QVector<ScoredCandidate> scoredCandidatesClassifier;
    e = fdrClassifierNeuralNet.exec(
            keyVsScoredCandidateCulled,
            scoredCandidatesAllFullFragIons,
            &scoredCandidatesClassifier
            );
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FDRCLassifierNeuralNetTests)
#include "FDRCLassifierNeuralNetTests.moc"
