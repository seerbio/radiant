#include "PeptideStringWithMods.h"

#include "AminoAcids.h"
#include "Molecule.h"

#include <QtTest>


class PeptideStringWithModsTests : public QObject
{
    Q_OBJECT

public:
    PeptideStringWithModsTests() = default;
    ~PeptideStringWithModsTests() override = default;

private slots:

    static void setTest();
    static void removeUniModCharsTest();
    static void modificationsMapTest();
    static void sizeNoModsTest();
    static void bSeriesTest();
    static void ySeriesTest();
    static void bSeriesIonLabelsTest();
    static void ySeriesIonLabelsTest();
    static void isValidSequenceTest();

};

void PeptideStringWithModsTests::setTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCYANDFLOPS");
    QCOMPARE(peptideStringWithMods, PeptideStringWithMods("CHAUNCYANDFLOPS"));
}

void PeptideStringWithModsTests::removeUniModCharsTest() {
    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCY[666]FLOPS");
    QCOMPARE(peptideStringWithMods.removeUniModChars(), QStringLiteral("CHAUNCYFLOPS"));
}

void PeptideStringWithModsTests::modificationsMapTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("B[66.6]ELL[+6666.6]ATRIX[-666.6]");
    const QMap<Index, double> modsMap = peptideStringWithMods.modificationsMap();

    QCOMPARE(modsMap.size(), 3);
    QCOMPARE(modsMap.firstKey(), 0);
    QCOMPARE(modsMap.first(), 66.6);
    QCOMPARE(modsMap.lastKey(), 8);
    QCOMPARE(modsMap.last(), -666.6);
    QCOMPARE(modsMap.keys().at(1), 3);
    QCOMPARE(modsMap.values().at(1), 6666.6);

    const PeptideStringWithMods peptideStringWithModsUniMod = PeptideStringWithMods("K(UniMod:4)ALL(Carbamidomethyl)IOPE(CAM)");
    const QMap<Index, double> modsMapUniMod = peptideStringWithModsUniMod.modificationsMap();

    QCOMPARE(modsMapUniMod.size(), 3);
    QCOMPARE(modsMapUniMod.firstKey(), 0);
    QCOMPARE(static_cast<int>(modsMapUniMod.first()), 57);
    QCOMPARE(modsMapUniMod.lastKey(), 7);
    QCOMPARE(static_cast<int>(modsMapUniMod.last()), 57);
    QCOMPARE(modsMapUniMod.keys().at(1), 3);
    QCOMPARE(static_cast<int>(modsMapUniMod.values().at(1)), 57);

    const PeptideStringWithMods peptideStringWithModsNTerm = PeptideStringWithMods("(UniMod:1)ACD");
    const QMap<Index, double> modsMapNTerm = peptideStringWithModsNTerm.modificationsMap();

    QCOMPARE(modsMapNTerm.size(), 1);
    QCOMPARE(modsMapNTerm.firstKey(), -1);
    QCOMPARE(MathUtils::pRound(modsMapNTerm.first(), 6), MathUtils::pRound(42.010565, 6));

    const PeptideStringWithMods peptideStringWithModsCTerm = PeptideStringWithMods("ACD(42.010565)");
    const QMap<Index, double> modsMapCTerm = peptideStringWithModsCTerm.modificationsMap();

    QCOMPARE(modsMapCTerm.size(), 1);
    QCOMPARE(modsMapCTerm.firstKey(), 2);
    QCOMPARE(MathUtils::pRound(modsMapCTerm.first(), 6), MathUtils::pRound(42.010565, 6));

}

void PeptideStringWithModsTests::sizeNoModsTest() {
    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCY[666]FLOPS");
    QCOMPARE(peptideStringWithMods.sizeNoMods(), 12);
}

void PeptideStringWithModsTests::bSeriesTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("V[1.0]T[-1.0]V[1.0]T[-1.0]");

    const QVector<double> bSeriesCharge1 = peptideStringWithMods.bSeries(1, AminoAcids());
    const QVector<double> expectedCharge1 = {101.076, 201.123, 301.192, 419.25};
    for (int i = 0; i < bSeriesCharge1.size(); i++) {
        QCOMPARE(MathUtils::pRound(bSeriesCharge1.at(i), 3), MathUtils::pRound(expectedCharge1.at(i), 3));
    }

    const QVector<double> bSeriesCharge2 = peptideStringWithMods.bSeries(2, AminoAcids());
    const QVector<double> expectedCharge2 = {51.041, 101.065, 151.1, 210.129};
    for (int i = 0; i < bSeriesCharge2.size(); i++) {
        QCOMPARE(MathUtils::pRound(bSeriesCharge2.at(i), 3), MathUtils::pRound(expectedCharge2.at(i), 3));
    }

}

void PeptideStringWithModsTests::ySeriesTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("V[1.0]T[-1.0]V[1.0]T[-1.0]");

    const QVector<double> ySeriesCharge1 = peptideStringWithMods.ySeries(1, AminoAcids());
    const QVector<double> expectedCharge1 = {119.066, 219.134, 319.182, 419.25};
    for (int i = 0; i < ySeriesCharge1.size(); i++) {
        QCOMPARE(MathUtils::pRound(ySeriesCharge1.at(i), 3), MathUtils::pRound(expectedCharge1.at(i), 3));
    }

    const QVector<double> ySeriesCharge2 = peptideStringWithMods.ySeries(2, AminoAcids());
    const QVector<double> expectedCharge2 = {60.0364, 110.071, 160.094, 210.129};
    for (int i = 0; i < ySeriesCharge2.size(); i++) {
        QCOMPARE(MathUtils::pRound(ySeriesCharge2.at(i), 3), MathUtils::pRound(expectedCharge2.at(i), 3));
    }

}

void PeptideStringWithModsTests::bSeriesIonLabelsTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("V[1.0]T[-1.0]V[1.0]T[-1.0]");
    const QStringList expected = {"b1", "b2", "b3", "b4"};
    QCOMPARE(peptideStringWithMods.bSeriesIonLabels({}), expected);

    const QStringList expected2 = {"b1^2", "b2^2", "b3^2", "b4^2"};
    QCOMPARE(peptideStringWithMods.bSeriesIonLabels("^2"), expected2);

}

void PeptideStringWithModsTests::ySeriesIonLabelsTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("V[1.0]T[-1.0]V[1.0]T[-1.0]");
    const QStringList expected = {"y1", "y2", "y3", "y4"};
    QCOMPARE(peptideStringWithMods.ySeriesIonLabels({}), expected);

    const QStringList expected2 = {"y1^2", "y2^2", "y3^2", "y4^2"};
    QCOMPARE(peptideStringWithMods.ySeriesIonLabels("^2"), expected2);

}

void PeptideStringWithModsTests::isValidSequenceTest() {
    const PeptideStringWithMods valid = PeptideStringWithMods("JOA[+15.99]CDEFGHIKLM(Oxidation)NPQRSTVWY");
    QCOMPARE(valid.isValidSequence(), true);

    PeptideStringWithMods invalid = PeptideStringWithMods("A[+15.99]PEPWITH(Oxidation)XK");
    QCOMPARE(invalid.isValidSequence(), false);
    invalid = PeptideStringWithMods("A[+15.99]PEPWITH(Oxidation)UK");
    QCOMPARE(invalid.isValidSequence(), false);
    invalid = PeptideStringWithMods("A[+15.99]PEPWITH(Oxidation)BK");
    QCOMPARE(invalid.isValidSequence(), false);
    invalid = PeptideStringWithMods("A[+15.99]PEPWITH(Oxidation)ZK");
    QCOMPARE(invalid.isValidSequence(), false);
}


QTEST_MAIN(PeptideStringWithModsTests)

#include "PeptideStringWithModsTests.moc"
