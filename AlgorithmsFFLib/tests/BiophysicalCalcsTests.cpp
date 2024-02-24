#include "BiophysicalCalcs.h"

#include "AminoAcids.h"
#include "PeptideStringWithMods.h"

#include <QHash>
#include <QtTest>

class BiophysicalCalcsTests : public QObject
{
    Q_OBJECT

public:
    BiophysicalCalcsTests() = default;
    ~BiophysicalCalcsTests() override = default;

private slots:

    static void calculatePeptideMassTest();
    static void calculateMassFromThomsonTest();
    static void calculatePeptideMassesTest();
    static void buildTandemFragmentMassesTest();
    static void calculateThomsonTest();
    static void calculateChargeFromSequenceTest();
    static void calculateIsotopesFromMzTest();

};

void BiophysicalCalcsTests::calculatePeptideMassTest() {

    const QString seq = QStringLiteral("TT");
    const QHash<int, double> mod = {
            {0, 100.0}
    };

    const double seqMass = BiophysicalCalcs::calculatePeptideMass(
            seq,
            AminoAcids(),
            mod
            );

    const double expectedSeqMass = 320.106;

    const auto test = static_cast<int>(seqMass);
    const auto expected = static_cast<int>(expectedSeqMass);

    QCOMPARE(test, expected);
}


void BiophysicalCalcsTests::calculatePeptideMassesTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("AMK");
    const QHash<ResidueIndex, ModificationMass> mods = {{'M', 16.0}};

    QVector<QPair<PeptideString, QHash<ResidueIndex, ModificationMass>>> vec = {
            {peptideStringWithMods, mods}
    };

    QVector<QPair<PeptideString, double>> result = BiophysicalCalcs::calculatePeptideMasses(
            vec,
            AminoAcids()
            );

    QCOMPARE(result.front().first, peptideStringWithMods);
    QCOMPARE(MathUtils::pRound(result.front().second, 3), 364.183);

}

void BiophysicalCalcsTests::calculateMassFromThomsonTest() {

    const double mass = BiophysicalCalcs::calculateMassFromThomson(1000.0, 2, 0);
    QCOMPARE(mass, 1997.98544706);

    const double massOffset = BiophysicalCalcs::calculateMassFromThomson(1000.0, 2, 1);
    QCOMPARE(massOffset, 1996.98209206);

}

void BiophysicalCalcsTests::buildTandemFragmentMassesTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("AMK");

    const QVector<double> frags = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideStringWithMods,
            2,
            100.0,
            CommonMolecules::H2O.monoisotopicMass(),
            AminoAcids()
            );


    const QVector<double> expected =  {86.5258, 152.046, 225.099};

    QCOMPARE(frags.size(), expected.size());

    for (int i = 0; i < frags.size(); i++) {
        QCOMPARE(MathUtils::pRound(frags.at(i), 3), MathUtils::pRound(expected.at(i), 3));
    }

    const PeptideStringWithMods peptideStringWithModsModsAdded = PeptideStringWithMods("V[+1.0]T[-1.0]V[+1.0]VVVV[-1.0]");

    const QVector<double> fragsWithMods = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideStringWithModsModsAdded,
            1,
            0.5,
            -0.5,
            AminoAcids()
    );
    
    const QVector<double> expectedModsAdded = {101.576, 201.623, 301.692, 400.76, 499.829, 598.897, 696.465};
    QCOMPARE(fragsWithMods.size(), expectedModsAdded.size());

    for (int i = 0; i < fragsWithMods.size(); i++) {
        QCOMPARE(MathUtils::pRound(fragsWithMods.at(i), 3), MathUtils::pRound(expectedModsAdded.at(i), 3));
    }

}

void BiophysicalCalcsTests::calculateThomsonTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("AMK");
    const double thomson = BiophysicalCalcs::calculateThomson(peptideStringWithMods, AminoAcids(), 2);
    QCOMPARE(thomson, 175.09883997);
}

void BiophysicalCalcsTests::calculateChargeFromSequenceTest() {

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("AMK");

    const int charge = BiophysicalCalcs::calculateChargeFromSequence(
            peptideStringWithMods,
            AminoAcids(),
            175.0988
            );
    QCOMPARE(charge, 2);

}

void BiophysicalCalcsTests::calculateIsotopesFromMzTest() {

    const QVector<double> isotopes = BiophysicalCalcs::calculateIsotopesFromMz(
            1000.0,
            2,
            1,
            2
            );

    const QVector<double> expected =  {999.498, 1000.000, 1000.502, 1001.003};

    QCOMPARE(isotopes.size(), expected.size());

    for (int i = 0; i < isotopes.size(); i++) {
        QCOMPARE(MathUtils::pRound(isotopes.at(i), 3), MathUtils::pRound(expected.at(i), 3));
    }

}


QTEST_MAIN(BiophysicalCalcsTests)

#include "BiophysicalCalcsTests.moc"
