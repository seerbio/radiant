#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "Molecule.h"
#include "MolecularFormula.h"

#include "MathUtils.h"

using namespace Error;

class MoleculeTests : public QObject
{
    Q_OBJECT

public:
    MoleculeTests()
    : m_testFormula(MolecularFormula(5, 9, 1, 1, 1, 0)){}
    ~MoleculeTests() override = default;

private slots:
    void averageMassTest();
    void monoisotopicMass();
    void parseMolecularFormulaStringTest();

private:
    MolecularFormula m_testFormula;

};

void MoleculeTests::averageMassTest() {
    Molecule testMolecule(m_testFormula);
    QCOMPARE(MathUtils::pRound(testMolecule.averageMass(), 3), 131.196);
}

void MoleculeTests::monoisotopicMass() {
    Molecule testMolecule(m_testFormula);
    QCOMPARE(MathUtils::pRound(testMolecule.monoisotopicMass(), 5), 131.04048);
}

void MoleculeTests::parseMolecularFormulaStringTest() {

    ERR_INIT

    const QString molString = QStringLiteral("H(-1)N-1O");

    MolecularFormula mf;

    e = MolecularFormulas::parseMolecularFormulaString(molString, &mf);
    QCOMPARE(e, eNoError);
    Molecule deamidation(mf);
    QCOMPARE(QString::number(deamidation.monoisotopicMass()), QString::number(0.984016));

    const QString molString2 = QStringLiteral("H(-3) C(-2) N(-1) O (-1) S");
    e = MolecularFormulas::parseMolecularFormulaString(molString2, &mf);
    QCOMPARE(e, eNoError);
    Molecule glnCysSubstitution(mf);
    QCOMPARE(QString::number(glnCysSubstitution.monoisotopicMass()), QString::number(-25.049393));

    const QString molString3 = QStringLiteral("H(-9) C(-5) N(-1) O(-1) S(-1)");
    e = MolecularFormulas::parseMolecularFormulaString(molString3, &mf);
    QCOMPARE(e, eNoError);
    Molecule metLoss(mf);
    QCOMPARE(QString::number(metLoss.monoisotopicMass()), QString::number(-131.040485));

    const QString errorMol = QStringLiteral("Q(-9) C(-5) N(-1) O(-1) S(-1)");
    e = MolecularFormulas::parseMolecularFormulaString(errorMol, &mf);
    QCOMPARE(e, eError);

}


QTEST_MAIN(MoleculeTests)

#include "MoleculeTests.moc"
