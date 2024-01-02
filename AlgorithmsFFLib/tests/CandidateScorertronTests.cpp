#include <QtTest>
#include <QCoreApplication>

#include "CandidateScorertron.h"

#include <QVector>

class CandidateScorertronTests : public QObject
{
    Q_OBJECT

public:
    CandidateScorertronTests() = default;
    ~CandidateScorertronTests() override = default;

private slots:
    void getIsotopicDistributionTest1();

};


void CandidateScorertronTests::getIsotopicDistributionTest1()
{




}



QTEST_MAIN(CandidateScorertronTests)

#include "CandidateScorertronTests.moc"
