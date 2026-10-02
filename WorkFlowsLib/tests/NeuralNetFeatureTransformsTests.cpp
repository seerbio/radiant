#include "NeuralNetFeatureTransforms.h"

#include <QtTest>

class NeuralNetFeatureTransformsTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void keepsOriginVariantsInOneFold();
    void transformsOnlyRawIntensities();
};

void NeuralNetFeatureTransformsTests::keepsOriginVariantsInOneFold() {
    const QStringList origins = {
        "PEPTIDEK", "_PEPTIDEK_", "PEPTLDEK", "_PEPT(Phospho)IDEK_"
    };
    const auto expected = NeuralNetFeatureTransforms::peptideFamilyHash(origins.first());
    for (const QString &origin : origins) {
        QCOMPARE(NeuralNetFeatureTransforms::peptideFamily(origin), QString("PEPTLDEK"));
        QCOMPARE(NeuralNetFeatureTransforms::peptideFamilyHash(origin), expected);
    }
    QVERIFY(expected != NeuralNetFeatureTransforms::peptideFamilyHash("ACDEFGHK"));
}

void NeuralNetFeatureTransformsTests::transformsOnlyRawIntensities() {
    const QVector<Features> features = {
        TotalIntensityRaw, Ms1IntensityFound100, Ms1IntensityFoundIso1,
        Ms1IntensityFoundPreMono, CosineSim100MS1, TotalIntensityLog, DiscriminantScore
    };
    QVector<float> values = {0.0f, 100.0f, 10000.0f, -1.0f, .75f, 7.0f, -3.0f};
    NeuralNetFeatureTransforms::logIntensities(features, &values);
    QCOMPARE(values[0], 0.0f);
    QVERIFY(std::abs(values[1] - std::log(101.0f)) < 1e-6f);
    QVERIFY(std::abs(values[2] - std::log(10001.0f)) < 1e-6f);
    QCOMPARE(values[3], 0.0f);
    QCOMPARE(values[4], .75f);
    QCOMPARE(values[5], 7.0f);
    QCOMPARE(values[6], -3.0f);
}

QTEST_MAIN(NeuralNetFeatureTransformsTests)
#include "NeuralNetFeatureTransformsTests.moc"
