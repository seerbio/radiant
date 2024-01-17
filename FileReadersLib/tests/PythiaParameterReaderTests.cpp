//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "PythiaParameterReader.h"

#include <QtTest/QtTest>
#include <QVariant>

class PythiaParameterReaderTests : public QObject
{
    Q_OBJECT

public:
    PythiaParameterReaderTests() = default;
    ~PythiaParameterReaderTests() override = default;

private Q_SLOTS:

    void readFileTest();

};

void PythiaParameterReaderTests::readFileTest() {

    ERR_INIT

    const QString filePath
            = QDir(qApp->applicationDirPath()).filePath("test_params_wide_window.pythiaConfig");

    PythiaParameterReader reader;

    PythiaParameters pythiaParameters;
    PythiaParameterReader::buildPythiaParameters(filePath, &pythiaParameters);
    QCOMPARE(e, eNoError);

    const QStringList cTermExpected = {"K", "R"};
    const QStringList nTermExpected = {};

    QCOMPARE(pythiaParameters.cTermCleavePoints, cTermExpected);
    QCOMPARE(pythiaParameters.nTermCleavePoints, nTermExpected);
    QCOMPARE(pythiaParameters.raggedness, 0);
    QCOMPARE(pythiaParameters.allowedMissedCleavages, 1);
    QCOMPARE(pythiaParameters.peptideLengthMin, 7);
    QCOMPARE(pythiaParameters.peptideLengthMax, 35);
    QCOMPARE(pythiaParameters.maxModificationsPeptide, 1);

    QCOMPARE(pythiaParameters.minFoundMzPeaks, 2);
    QCOMPARE(pythiaParameters.ms2ExtractionWidthPPM, 20.5);
    QCOMPARE(pythiaParameters.mzMinMS2, 100.1);
    QCOMPARE(pythiaParameters.mzMaxMS2, 1500.2);

    QCOMPARE(pythiaParameters.chargeStateMin, 1);
    QCOMPARE(pythiaParameters.chargeStateMax, 4);
    QCOMPARE(pythiaParameters.precursorExtractionWindowThomsons, 0.5);

    QCOMPARE(pythiaParameters.filterLength, 3);
    QCOMPARE(pythiaParameters.sigma, 1.0);
    QCOMPARE(pythiaParameters.signalToNoiseRatio, 2.0);
    QCOMPARE(pythiaParameters.smoothCount, 1);

    QCOMPARE(pythiaParameters.percentFDR, 2.0);
    QCOMPARE(pythiaParameters.reportDecoys, true);

    QCOMPARE(pythiaParameters.modifications.size(), 3);

    const Modification &modFixed = pythiaParameters.modifications.at(0);
    QCOMPARE(modFixed.name, "Carbamidomethyl");
    QCOMPARE(modFixed.type, ModificationType::FIXED);
    QCOMPARE(modFixed.formula, "H3C2NO");
    QCOMPARE(modFixed.residue, "C");
    QCOMPARE(
            QString::number(pythiaParameters.aminoAcids.aminoAcid('C').monoisotopicMass()),
            QString::number(160.031)
            );

    const Modification &modVar1 = pythiaParameters.modifications.at(1);
    QCOMPARE(modVar1.name, "Oxidation");
    QCOMPARE(modVar1.type, ModificationType::DYNAMIC);
    QCOMPARE(modVar1.formula, "O");
    QCOMPARE(modVar1.residue, "M");

    MolecularFormula mfOx;
    e = parseMolecularFormulaString(modVar1.formula, &mfOx);
    QCOMPARE(e, eNoError);
    QCOMPARE(QString::number(Molecule(mfOx).monoisotopicMass()), QString::number(15.9949));

    const Modification &modVar2 = pythiaParameters.modifications.at(2);
    QCOMPARE(modVar2.name, "Acetyl");
    QCOMPARE(modVar2.type, ModificationType::DYNAMIC);
    QCOMPARE(modVar2.formula, "C2H3O");
    QCOMPARE(modVar2.positionalLocation, "N-Term-Protein");

    MolecularFormula mfAcetyl;
    e = parseMolecularFormulaString(modVar2.formula, &mfAcetyl);
    QCOMPARE(e, eNoError);
    QCOMPARE(QString::number(Molecule(mfAcetyl).monoisotopicMass()), QString::number(43.0184));
}


QTEST_MAIN(PythiaParameterReaderTests)
#include "PythiaParameterReaderTests.moc"
