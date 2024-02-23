#include "PeptideStringWithMods.h"
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

};

void PeptideStringWithModsTests::setTest() {

    PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCYANDFLOPS");
    QCOMPARE(peptideStringWithMods, PeptideStringWithMods("CHAUNCYANDFLOPS"));
}

void PeptideStringWithModsTests::removeUniModCharsTest() {
    PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCY[666]FLOPS");
    QCOMPARE(peptideStringWithMods.removeUniModChars(), QStringLiteral("CHAUNCYFLOPS"));
}

void PeptideStringWithModsTests::modificationsMapTest() {

    PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("B[66.6]ELL[+6666.6]ATRIX[-666.6]");
    const QMap<Index, Mass> modsMap = peptideStringWithMods.modificationsMap();

    QCOMPARE(modsMap.size(), 3);
    QCOMPARE(modsMap.firstKey(), 0);
    QCOMPARE(modsMap.first(), 66.6);
    QCOMPARE(modsMap.lastKey(), 8);
    QCOMPARE(modsMap.last(), -666.6);
    QCOMPARE(modsMap.keys().at(1), 3);
    QCOMPARE(modsMap.values().at(1), 6666.6);

    PeptideStringWithMods peptideStringWithModsUniMod = PeptideStringWithMods("K(UniMod:4)ALL(Carbamidomethyl)IOPE(CAM)");
    const QMap<Index, Mass> modsMapUniMod = peptideStringWithModsUniMod.modificationsMap();

    QCOMPARE(modsMapUniMod.size(), 3);
    QCOMPARE(modsMapUniMod.firstKey(), 0);
    QCOMPARE(static_cast<int>(modsMapUniMod.first()), 57);
    QCOMPARE(modsMapUniMod.lastKey(), 7);
    QCOMPARE(static_cast<int>(modsMapUniMod.last()), 57);
    QCOMPARE(modsMapUniMod.keys().at(1), 3);
    QCOMPARE(static_cast<int>(modsMapUniMod.values().at(1)), 57);


}

void PeptideStringWithModsTests::sizeNoModsTest() {
    PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCY[666]FLOPS");
    QCOMPARE(peptideStringWithMods.sizeNoMods(), 12);
}


QTEST_MAIN(PeptideStringWithModsTests)

#include "PeptideStringWithModsTests.moc"
