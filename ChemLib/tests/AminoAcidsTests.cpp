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
}

QTEST_MAIN(AminoAcidsTests)

#include "AminoAcidsTests.moc"
