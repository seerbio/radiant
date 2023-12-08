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

    const PeptideStringWithMods peptideStringWithMods = QStringLiteral("CHAUNCYANDFLOPS");
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
            totalFramentCount
            );

    QCOMPARE(targetDecoyCandidatePair.peptideStringWithMods(), peptideStringWithMods);
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsTarget().first().mz, 666.6));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsDecoy().first().mz, 777.7));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mz(), 334.3072));
    QCOMPARE(targetDecoyCandidatePair.charge(), charge);
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mass(), 666.6));
    QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.iRt(), 66.6));
    QCOMPARE(targetDecoyCandidatePair.totalFragmentCount(), 12);

}


QTEST_MAIN(TargetDecoyCandidatePairTests)
#include "TargetDecoyCandidatePairTests.moc"
