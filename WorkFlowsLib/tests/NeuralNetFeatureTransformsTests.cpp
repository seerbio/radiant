#include "NeuralNetFeatureTransforms.h"
#include "DiscriminantScoretron.h"

#include <QtTest>

#include <cmath>

class NeuralNetFeatureTransformsTests : public QObject {
Q_OBJECT

private slots:

    static void logsOnlyRawMs1IntensityFeatures();
    static void neuralNetworkFeaturesUseLogTotalIntensity();
};

void NeuralNetFeatureTransformsTests::logsOnlyRawMs1IntensityFeatures() {

    const QVector<Features> features = {
        TotalIntensityLog,
        TotalIntensityRaw,
        Ms1IntensityFound100,
        Ms1IntensityFound45,
        Ms1IntensityFoundPreMono,
        Ms1IntensityFoundIso1,
        Ms1IntensityFoundIso2,
        Ms1IntensityFoundApex100IM,
        CosineSim100MS1
    };
    QVector<float> values = {
        7.0f,
        1000.0f,
        0.0f,
        100.0f,
        10000.0f,
        -1.0f,
        25.0f,
        50.0f,
        0.75f
    };

    NeuralNetFeatureTransforms::logIntensityFeatures(features, &values);

    QCOMPARE(values.at(0), 7.0f);
    QCOMPARE(values.at(1), 1000.0f);
    QVERIFY(std::abs(values.at(2) - std::log1p(0.0f)) < 1e-6f);
    QVERIFY(std::abs(values.at(3) - std::log1p(100.0f)) < 1e-6f);
    QVERIFY(std::abs(values.at(4) - std::log1p(10000.0f)) < 1e-6f);
    QVERIFY(std::abs(values.at(5) - std::log1p(0.0f)) < 1e-6f);
    QVERIFY(std::abs(values.at(6) - std::log1p(25.0f)) < 1e-6f);
    QVERIFY(std::abs(values.at(7) - std::log1p(50.0f)) < 1e-6f);
    QCOMPARE(values.at(8), 0.75f);
}

void NeuralNetFeatureTransformsTests::neuralNetworkFeaturesUseLogTotalIntensity() {

    const QVector<Features> features = DiscriminantScoretron::featuresNeuralNetwork();
    QVERIFY(features.contains(TotalIntensityLog));
    QVERIFY(!features.contains(TotalIntensityRaw));
}

QTEST_MAIN(NeuralNetFeatureTransformsTests)

#include "NeuralNetFeatureTransformsTests.moc"
