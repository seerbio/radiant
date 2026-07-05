#include <QtTest>
#include <QCoreApplication>

#include "CandidateScores.h"
#include "DiscriminantScoretron.h"

class DiscriminantScoretronTests : public QObject {
Q_OBJECT

public:
    DiscriminantScoretronTests() = default;

    ~DiscriminantScoretronTests() override = default;

private slots:

    static void setDisciminantScoreForCandidatesTest();
    static void defaultWeightsPreserveLegacyNonIonMobilityFeaturesTest();
    static void defaultWeightsEnableIonMobilityFeaturesWhenSelectedTest();



};

void DiscriminantScoretronTests::setDisciminantScoreForCandidatesTest() {
    QSKIP("Legacy no-op fixture; covered by default weight tests.");
}

void DiscriminantScoretronTests::defaultWeightsPreserveLegacyNonIonMobilityFeaturesTest() {

    const QVector<Features> features = {
        CosineSimSum100GreaterThan80,
        CosineSimSum100Top12,
        CosineSimSpectrumOverTimeCubed,
        ScanTimeDeltaAbs
    };

    const QVector<float> weights = DiscriminantScoretron::defaultWeights(features);

    QCOMPARE(weights.size(), features.size());
    QCOMPARE(weights.at(0), 1.0f);
    QCOMPARE(weights.at(1), 1.0f);
    QCOMPARE(weights.at(2), 1.0f);
    QCOMPARE(weights.at(3), -0.5f);
}

void DiscriminantScoretronTests::defaultWeightsEnableIonMobilityFeaturesWhenSelectedTest() {

    const QVector<Features> features = {
        IonMobilityDeltaAbs,
        IonMobilityPdAbs,
        Ms2IonMobilityWeightedDeltaAbs,
        Ms2IonMobilityMatchedIonFraction
    };

    const QVector<float> weights = DiscriminantScoretron::defaultWeights(features);

    QCOMPARE(weights.size(), features.size());
    QCOMPARE(weights.at(0), -1.0f);
    QCOMPARE(weights.at(1), -0.5f);
    QCOMPARE(weights.at(2), -1.0f);
    QCOMPARE(weights.at(3), 1.0f);
}


QTEST_MAIN(DiscriminantScoretronTests)

#include "DiscriminantScoretronTests.moc"
