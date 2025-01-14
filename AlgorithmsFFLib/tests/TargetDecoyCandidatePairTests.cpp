#include "TargetDecoyCandidatePair.h"

#include <QtTest/QtTest>


class TargetDecoyCandidatePairTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairTests() = default;
    ~TargetDecoyCandidatePairTests() override = default;

private Q_SLOTS:

    void gettersTest();


};

void TargetDecoyCandidatePairTests::gettersTest() {

    MS2Ion ms2Ion1;
    ms2Ion1.mz = 666.6;

    MS2Ion ms2Ion2;
    ms2Ion2.mz = 777.7;

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCYANDFLOPS");
    const QVector<MS2Ion> ms2IonsTarget = {ms2Ion1};
    const QVector<MS2Ion> ms2IonsDecoy = {ms2Ion2};
    const int charge = 2;
    const double mass = 666.6;
    const double iRt = 66.6;
    const int totalFramentCount = 12;

    TargetDecoyCandidatePair targetDecoyCandidatePair(
            peptideStringWithMods,
            ms2IonsTarget,
            ms2IonsDecoy,
            charge,
            mass,
            iRt,
            iRt,
            totalFramentCount,
            0.0
            );

    QCOMPARE(targetDecoyCandidatePair.peptideStringWithMods(), peptideStringWithMods);
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsTarget().first().mz, 666.6f));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsDecoy().first().mz, 777.7f));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mz(false), 334.3072f));
    QCOMPARE(targetDecoyCandidatePair.charge(), charge);
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mass(), 666.6f));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.iRt(), 66.6f));
    QCOMPARE(targetDecoyCandidatePair.totalFragmentCount(), 12);

}


QTEST_MAIN(TargetDecoyCandidatePairTests)
#include "TargetDecoyCandidatePairTests.moc"
