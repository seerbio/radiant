#include <QtTest>
#include <QCoreApplication>

#include "IsotopicDistributionBuilder.h"
#include "MathUtils.h"

class IsotopicDistributionBuilderTests : public QObject {
Q_OBJECT

public:
    IsotopicDistributionBuilderTests() = default;

    ~IsotopicDistributionBuilderTests() override = default;

private slots:

    static void getIsotopicDistributionTest1();

    static void getIsotopicDistributionTest2();

    static void getMolecularFormulaOfBioPolymerTest();



};

void IsotopicDistributionBuilderTests::getIsotopicDistributionTest1() {

    ERR_INIT

    const QString seq = QStringLiteral("IHEARTCATSIHEARTCATS");
    const QStringList seqList = seq.split("");

    const QVector<double> isoDist = IsotopicDistributionBuilder::getIsotopicDistribution(seqList);
    const QVector<double> expectedIsoDist = {0.917159, 1, 0.639969, 0.301553, 0.113406, 0.0355548, 0.0095633, 0.0022531};

    QCOMPARE(isoDist.size(), expectedIsoDist.size());
    for (int i = 0; i < isoDist.size(); i++) {
        QCOMPARE(MathUtils::pRound(isoDist.at(i), 3), MathUtils::pRound(expectedIsoDist.at(i), 3));
    }

}

void IsotopicDistributionBuilderTests::getIsotopicDistributionTest2() {

    ERR_INIT

    const QVector<double> isoDist = IsotopicDistributionBuilder::getIsotopicDistribution(1800.0);
    const QVector<double> expectedIsoDist = {1, 0.990332, 0.548216, 0.219205, 0.0691112, 0.0179445, 0.00395654, 0.000760633};

    QCOMPARE(isoDist.size(), expectedIsoDist.size());
    for (int i = 0; i < isoDist.size(); i++) {
        QCOMPARE(MathUtils::pRound(isoDist.at(i), 3), MathUtils::pRound(expectedIsoDist.at(i), 3));
    }

}

void IsotopicDistributionBuilderTests::getMolecularFormulaOfBioPolymerTest() {

    ERR_INIT

    const QString seq = QStringLiteral("IHEARTCATSIHEARTCATS");
    const QStringList seqList = seq.split("");

    const MolecularFormula mf = IsotopicDistributionBuilder::getMolecularFormulaOfBioPolymer(seqList);
    QCOMPARE(mf.carbonCount, 86);
    QCOMPARE(mf.hydrogenCount, 144);
    QCOMPARE(mf.nitrogenCount, 30);
    QCOMPARE(mf.oxygenCount, 31);
    QCOMPARE(mf.sulfurCount, 2);
    QCOMPARE(mf.phosphorusCount, 0);

}


QTEST_MAIN(IsotopicDistributionBuilderTests)

#include "IsotopicDistributionBuilderTests.moc"
