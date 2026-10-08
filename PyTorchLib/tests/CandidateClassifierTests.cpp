//
// Created by anichols on 11/07/2021.
//

#include "CandidateClassifier.h"

#include <QDebug>
#include <QtTest/QtTest>


class CandidateClassifierTests : public QObject
{
    Q_OBJECT

public:
    CandidateClassifierTests() = default;
    ~CandidateClassifierTests() override = default;

private Q_SLOTS:

    static void trainCandidateClassifierAndPredictTest();
    static void predictionDoesNotChangeSubsequentTraining();

};

void CandidateClassifierTests::trainCandidateClassifierAndPredictTest() {

    const QVector<QVector<float>> xVec = {
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 1.0, 0.0}
    };

    const QVector<float> yVec = {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0};

    CandidateClassifier classifier;
    const bool classifierTrainedNoErrors = classifier.trainCandidateClassifier(xVec, yVec, 30, 2, 1e-2, 666, 0.5, 0.0, 0);
    QCOMPARE(classifierTrainedNoErrors, true);

    QVector<float> predictions;
    const bool predictedNoErrors = classifier.predict(xVec, &predictions);
    QCOMPARE(predictions.size(), yVec.size());
    QCOMPARE(predictedNoErrors, true);
    for (int i = 0; i < predictions.size(); i++) {
        qDebug() << predictions.at(i) << static_cast<int>(yVec.at(i));
        QCOMPARE(static_cast<int>(std::round(predictions.at(i))), static_cast<int>(yVec.at(i)));
    }

}

void CandidateClassifierTests::predictionDoesNotChangeSubsequentTraining() {
    const QVector<QVector<float>> xData = {
        {1.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 1.0f}, {0.0f, 0.0f},
        {1.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 1.0f}, {0.0f, 0.0f}
    };
    const QVector<float> yData = {1.0f, 0.0f, 1.0f, 0.0f,
                                  1.0f, 0.0f, 1.0f, 0.0f};
    const QVector<QVector<float>> inference = {
        {1.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 1.0f}
    };

    CandidateClassifier classifier;
    QVERIFY(classifier.trainCandidateClassifier(
        xData, yData, 10, 4, 0.003, 29, 1.0, 0.0f, 0));
    QVector<float> first, repeated, retrained;
    QVERIFY(classifier.predict(inference, &first));
    QVERIFY(classifier.predict(inference, &repeated));
    QCOMPARE(first, repeated);
    QVERIFY(classifier.trainCandidateClassifier(
        xData, yData, 10, 4, 0.003, 29, 1.0, 0.0f, 0));
    QVERIFY(classifier.predict(inference, &retrained));
    QCOMPARE(first, retrained);
}


QTEST_MAIN(CandidateClassifierTests)
#include "CandidateClassifierTests.moc"
