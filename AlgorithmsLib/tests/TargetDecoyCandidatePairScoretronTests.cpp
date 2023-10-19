#include "TargetDecoyCandidatePairScoretron.h"


#include <QtTest/QtTest>

#include <iostream>

class TargetDecoyCandidatePairScoretronTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairScoretronTests() = default;
    ~TargetDecoyCandidatePairScoretronTests() override = default;

private Q_SLOTS:
    void loadModelTest();


};

void TargetDecoyCandidatePairScoretronTests::loadModelTest() {

    ERR_INIT

    TargetDecoyCandidatePairScoretron targetDecoyCandidatePairScoretron;
//    e = targetDecoyCandidatePairScoretron.init();
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(TargetDecoyCandidatePairScoretronTests)
#include "TargetDecoyCandidatePairScoretronTests.moc"
