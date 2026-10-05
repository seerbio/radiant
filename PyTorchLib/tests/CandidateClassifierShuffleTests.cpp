#include "CandidateClassifier.h"

#include <QtTest>

#include <cmath>

namespace {

QVector<QVector<float>> observations(int count, float phase) {
    QVector<QVector<float>> result;
    result.reserve(count);
    for (int row = 0; row < count; ++row) {
        QVector<float> values(4);
        values[0] = row < count / 2 ? -0.8f : 0.8f;
        for (int column = 1; column < values.size(); ++column) {
            values[column] = 0.15f * std::sin(float(row * (column + 1)) + phase);
        }
        result.push_back(values);
    }
    return result;
}

QVector<float> labels(int count) {
    QVector<float> result(count);
    for (int row = 0; row < count; ++row) {
        result[row] = row < count / 2 ? 0.0f : 1.0f;
    }
    return result;
}

}

class CandidateClassifierShuffleTests : public QObject {
    Q_OBJECT

private slots:

    static void deterministicTrainingPreservesLabels();
    static void omittedOptionMatchesDisabled();
};

void CandidateClassifierShuffleTests::deterministicTrainingPreservesLabels() {
    const auto training = observations(600, 0.1f);
    const auto trainingLabels = labels(training.size());
    const auto inference = observations(120, 0.7f);
    const auto expected = labels(inference.size());

    QVector<QVector<float>> repeated;
    for (int repetition = 0; repetition < 2; ++repetition) {
        CandidateClassifier classifier;
        QVERIFY(classifier.trainCandidateClassifier(
            training, trainingLabels, 40, 50, 0.003, 17, 2.0, 0.0f, 0, true));

        QVector<float> probabilities;
        QVERIFY(classifier.predict(inference, &probabilities));
        QCOMPARE(probabilities.size(), expected.size());

        int correct = 0;
        for (int row = 0; row < probabilities.size(); ++row) {
            QVERIFY(std::isfinite(probabilities[row]));
            QVERIFY(probabilities[row] >= 0 && probabilities[row] <= 1);
            correct += (probabilities[row] >= 0.5f) == (expected[row] >= 0.5f);
        }
        QVERIFY2(correct >= 114, "Shuffling must preserve feature/label correspondence.");
        repeated.push_back(probabilities);
    }

    QCOMPARE(repeated[0], repeated[1]);
}

void CandidateClassifierShuffleTests::omittedOptionMatchesDisabled() {
    const auto training = observations(60, 0.2f);
    const auto trainingLabels = labels(training.size());
    CandidateClassifier original;
    CandidateClassifier explicitlyDisabled;

    QVERIFY(original.trainCandidateClassifier(
        training, trainingLabels, 4, 10, 0.003, 29, 1.0, 0.0f, 0));
    QVERIFY(explicitlyDisabled.trainCandidateClassifier(
        training, trainingLabels, 4, 10, 0.003, 29, 1.0, 0.0f, 0, false));

    QVector<float> originalPredictions;
    QVector<float> disabledPredictions;
    QVERIFY(original.predict(training, &originalPredictions));
    QVERIFY(explicitlyDisabled.predict(training, &disabledPredictions));
    QCOMPARE(originalPredictions, disabledPredictions);
}

QTEST_MAIN(CandidateClassifierShuffleTests)

#include "CandidateClassifierShuffleTests.moc"
