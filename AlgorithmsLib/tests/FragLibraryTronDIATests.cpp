#include "FragLibraryTronDIA.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "PythiaParameterReader.h"

#include <QtTest>

#include <iostream>

class FragLibraryTronDIATests : public QObject
{
    Q_OBJECT

public:

    FragLibraryTronDIATests() = default;
    ~FragLibraryTronDIATests() override = default;

private slots:

    void readTest();


private:


private:

    QString fragFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta.csv.fragLib");
    }

    PythiaParameters pythiaParameters() {

        PythiaParameters parms;
        parms.returnPSMTopN = 1;
        parms.maxTandemPointCount = 2;
        parms.ms2ExtractionWidthPPM = 12.0;
        parms.precursorExtractionWindowThomsons = 1.0;
        parms.chargeStateMin = 2;
        parms.chargeStateMax = 3;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                parms,
                &parms.aminoAcids
        );

        return parms;
    }

};

void FragLibraryTronDIATests::readTest()
{

    ERR_INIT;

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
            );
    QCOMPARE(e, eNoError);





}


QTEST_MAIN(FragLibraryTronDIATests)

#include "FragLibraryTronDIATests.moc"
