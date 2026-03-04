#include "AminoAcids.h"
#include "Molecule.h"

#include <QtTest>


class AminoAcidsTests : public QObject
{
    Q_OBJECT

public:
    AminoAcidsTests() = default;
    ~AminoAcidsTests() override = default;

private slots:
    void checkMassesTest();
    void mutatePenultimatePeptideResiduesTest();
    void validPeptideSequenceTest();

};


void AminoAcidsTests::checkMassesTest()
{

    const QMap<QChar, Molecule> AMINO_ACIDS = AminoAcids::aminoAcids();

    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('A').monoisotopicMass()), 71);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('C').monoisotopicMass()), 103);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('D').monoisotopicMass()), 115);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('E').monoisotopicMass()), 129);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('F').monoisotopicMass()), 147);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('G').monoisotopicMass()), 57);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('H').monoisotopicMass()), 137);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('I').monoisotopicMass()), 113);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('K').monoisotopicMass()), 128);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('L').monoisotopicMass()), 113);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('M').monoisotopicMass()), 131);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('N').monoisotopicMass()), 114);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('P').monoisotopicMass()), 97);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('Q').monoisotopicMass()), 128);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('R').monoisotopicMass()), 156);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('S').monoisotopicMass()), 87);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('T').monoisotopicMass()), 101);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('V').monoisotopicMass()), 99);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('W').monoisotopicMass()), 186);
    QCOMPARE(static_cast<int>(AMINO_ACIDS.value('Y').monoisotopicMass()), 163);

    QHash<QChar, double> AA;
    AA['G'] = 57.021464;
    AA['A'] = 71.037114;
    AA['V'] = 99.068414;
    AA['L'] = 113.084064;
    AA['I'] = 113.084064;

    AA['F'] = 147.068414;
    AA['M'] = 131.040485;
    AA['P'] = 97.052764;
    AA['W'] = 186.079313;
    AA['S'] = 87.032028;

    AA['C'] = 103.009185;
    AA['T'] = 101.047679;
    AA['Y'] = 163.06332;
    AA['H'] = 137.058912;
    AA['K'] = 128.094963;

    AA['R'] = 156.101111;
    AA['Q'] = 128.058578;
    AA['E'] = 129.042593;
    AA['N'] = 114.042927;
    AA['D'] = 115.026943;

//    qDebug() << AMINO_ACIDS.value('A').monoisotopicMass() - AA['A'];
//    qDebug() << AMINO_ACIDS.value('C').monoisotopicMass() - AA['C'];
//    qDebug() << AMINO_ACIDS.value('D').monoisotopicMass() - AA['D'];
//    qDebug() << AMINO_ACIDS.value('E').monoisotopicMass() - AA['E'];
//    qDebug() << AMINO_ACIDS.value('F').monoisotopicMass() - AA['F'];
//    qDebug() << AMINO_ACIDS.value('G').monoisotopicMass() - AA['G'];
//    qDebug() << AMINO_ACIDS.value('H').monoisotopicMass() - AA['H'];
//    qDebug() << AMINO_ACIDS.value('I').monoisotopicMass() - AA['I'];
//    qDebug() << AMINO_ACIDS.value('K').monoisotopicMass() - AA['K'];
//    qDebug() << AMINO_ACIDS.value('L').monoisotopicMass() - AA['L'];
//    qDebug() << AMINO_ACIDS.value('M').monoisotopicMass() - AA['M'];
//    qDebug() << AMINO_ACIDS.value('N').monoisotopicMass() - AA['N'];
//    qDebug() << AMINO_ACIDS.value('P').monoisotopicMass() - AA['P'];
//    qDebug() << AMINO_ACIDS.value('Q').monoisotopicMass() - AA['Q'];
//    qDebug() << AMINO_ACIDS.value('R').monoisotopicMass() - AA['R'];
//    qDebug() << AMINO_ACIDS.value('S').monoisotopicMass() - AA['S'];
//    qDebug() << AMINO_ACIDS.value('T').monoisotopicMass() - AA['T'];
//    qDebug() << AMINO_ACIDS.value('V').monoisotopicMass() - AA['V'];
//    qDebug() << AMINO_ACIDS.value('W').monoisotopicMass() - AA['W'];
//    qDebug() << AMINO_ACIDS.value('Y').monoisotopicMass() - AA['Y'];

}

void AminoAcidsTests::mutatePenultimatePeptideResiduesTest() {
    PeptideStringWithMods pep("SD(123)BC(Mod)D");

    PeptideStringWithMods modded = AminoAcids::mutatePenultimatePeptideResidues(pep);
    QCOMPARE(modded.at(1), "E");
    QCOMPARE(modded.at(8), "S");
}

void AminoAcidsTests::validPeptideSequenceTest() {
    QCOMPARE(AminoAcids::validPeptideSequence("ACDEFGHIKLMNPQRSTVWYJO"), true);
    QCOMPARE(AminoAcids::validPeptideSequence(""), true);

    QCOMPARE(AminoAcids::validPeptideSequence("B"), false);
    QCOMPARE(AminoAcids::validPeptideSequence("U"), false);
    QCOMPARE(AminoAcids::validPeptideSequence("X"), false);
    QCOMPARE(AminoAcids::validPeptideSequence("Z"), false);
    QCOMPARE(AminoAcids::validPeptideSequence("AJOX"), false);
}

QTEST_MAIN(AminoAcidsTests)

#include "AminoAcidsTests.moc"
