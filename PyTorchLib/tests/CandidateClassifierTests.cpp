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

    void readLibrary();

};

void CandidateClassifierTests::readLibrary() {

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
//    classifier.trainCandidateClassifier(xVec, yVec, 1000, 0.2, 1e-2);


}


QTEST_MAIN(CandidateClassifierTests)
#include "CandidateClassifierTests.moc"
