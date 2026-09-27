#include "CandidateClassifier.h"

#include <QtTest/QtTest>
#include <cmath>

class CandidateClassifierShuffleTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    static void deterministicTrainingPreservesLabels();
    static void omittedOptionMatchesDisabled();
};

namespace {
QVector<QVector<float>> observations(int count, float phase) {
    QVector<QVector<float>> result;
    for (int row = 0; row < count; ++row) {
        QVector<float> values(8);
        // Ordered classes make independent feature/label permutations fail.
        values[0] = row < count / 2 ? -.8f : .8f;
        for (int col = 1; col < values.size(); ++col)
            values[col] = .15f * std::sin(float(row * (col + 1)) + phase);
        result.push_back(values);
    }
    return result;
}
QVector<float> labels(int count) {
    QVector<float> result(count);
    for (int row = 0; row < count; ++row) result[row] = row < count / 2 ? 0.0f : 1.0f;
    return result;
}
}

void CandidateClassifierShuffleTests::deterministicTrainingPreservesLabels() {
    const auto training = observations(600, .1f);
    const auto y = labels(training.size());
    const auto inference = observations(120, .7f);
    const auto expected = labels(inference.size());
    QVector<QVector<float>> repeated;
    for (int repetition = 0; repetition < 2; ++repetition) {
        CandidateClassifier classifier;
        QVERIFY(classifier.trainCandidateClassifier(
            training, y, 40, 50, .003, 17, 2.0, 0.0f, 0, true));
        QVector<float> probabilities;
        QVERIFY(classifier.predict(inference, &probabilities));
        QCOMPARE(probabilities.size(), expected.size());
        int correct = 0;
        for (int row = 0; row < probabilities.size(); ++row) {
            QVERIFY(std::isfinite(probabilities[row]));
            QVERIFY(probabilities[row] >= 0 && probabilities[row] <= 1);
            correct += (probabilities[row] >= .5f) == (expected[row] >= .5f);
        }
        QVERIFY2(correct >= 114, "Shuffling must preserve feature/label correspondence.");
        repeated.push_back(probabilities);
    }
    QCOMPARE(repeated[0], repeated[1]);
}

void CandidateClassifierShuffleTests::omittedOptionMatchesDisabled() {
    const auto training = observations(60, .2f);
    const auto y = labels(training.size());
    CandidateClassifier original, explicitlyDisabled;
    QVERIFY(original.trainCandidateClassifier(training, y, 4, 10, .003, 29, 1.0, 0.0f, 0));
    QVERIFY(explicitlyDisabled.trainCandidateClassifier(training, y, 4, 10, .003, 29, 1.0, 0.0f, 0, false));
    QVector<float> a, b;
    QVERIFY(original.predict(training, &a));
    QVERIFY(explicitlyDisabled.predict(training, &b));
    QCOMPARE(a, b);
}

QTEST_MAIN(CandidateClassifierShuffleTests)
#include "CandidateClassifierShuffleTests.moc"
