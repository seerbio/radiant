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

    const QString &jsonFilePath
            = QDir(qApp->applicationDirPath()).filePath("DigestParamsTest2.pythia");


    PythiaParameterReader reader;
    e = reader.readFile(jsonFilePath);
    QCOMPARE(e, eNoError);

    PythiaParameters pythiaParameters;
    e = reader.loadPythiaParameters(&pythiaParameters);
    QCOMPARE(e, eNoError);

    pythiaParameters.print();

    QCOMPARE(pythiaParameters.nTermCleavePoints, QStringList({"S", "H", "I", "T"}));
    QCOMPARE(pythiaParameters.cTermCleavePoints, QStringList({"C", "R", "A", "P"}));

    QCOMPARE(pythiaParameters.raggedness, 3);
    QCOMPARE(pythiaParameters.allowedMissedCleavages, 22);
    QCOMPARE(pythiaParameters.mzMinDataStructure, 331);
    QCOMPARE(pythiaParameters.mzMaxDataStructure, 1992);
    QCOMPARE(pythiaParameters.peptideLengthMin, 8);
    QCOMPARE(pythiaParameters.peptideLengthMax, 30);
    QCOMPARE(pythiaParameters.maxTandemPointCount,401 );
    QCOMPARE(pythiaParameters.returnPSMTopN, 31);
    QCOMPARE(pythiaParameters.chargeStateMin, 5);
    QCOMPARE(pythiaParameters.chargeStateMax, 3);
    QCOMPARE(pythiaParameters.ms2ExtractionWidthPPM, 11);
    QCOMPARE(pythiaParameters.precursorExtractionWindowThomsons, 1);
    QCOMPARE(pythiaParameters.percentFDR, 2);
    QCOMPARE(pythiaParameters.maxModificationsPeptide, 5);
    QCOMPARE(pythiaParameters.addDecoys, true);

    const Modification &mod = pythiaParameters.modifications.at(0);
    QCOMPARE(mod.name, "Carbamidomethyl");
    QCOMPARE(mod.type, ModificationType::FIXED);
    QCOMPARE(mod.formula, "H3C2NO");
    QCOMPARE(mod.residue, "C");
}


QTEST_MAIN(PythiaParameterReaderTests)
#include "PythiaParameterReaderTests.moc"
