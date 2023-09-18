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

    CandidateClassifier r;


}


QTEST_MAIN(CandidateClassifierTests)
#include "CandidateClassifierTests.moc"
