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
    const bool classifierTrainedNoErrors = classifier.trainCandidateClassifier(xVec, yVec, 10, 2, 1e-2, 666);
    QCOMPARE(classifierTrainedNoErrors, true);

    QVector<float> predictions;
    const bool predictedNoErrors = classifier.predict(xVec, &predictions);
    QCOMPARE(predictions.size(), yVec.size());
    QCOMPARE(predictedNoErrors, true);
    for (int i = 0; i < predictions.size(); i++) {
        QCOMPARE(static_cast<int>(std::round(predictions.at(i))), static_cast<int>(yVec.at(i)));
    }

}


QTEST_MAIN(CandidateClassifierTests)
#include "CandidateClassifierTests.moc"
