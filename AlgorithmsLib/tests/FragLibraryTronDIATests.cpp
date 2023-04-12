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

    void initTest();

    void getMS2IonsTest();
    void getMS2IonsTest2();
    void getMS2IonsTest3();
    void getMS2IonsTest4();

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

void FragLibraryTronDIATests::initTest() {

//    QSKIP("TEMP");

    ERR_INIT;

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(fragLibraryTronDia.m_pepSeqChrgKeyVsMS2Ions.size(), 102362);

}

void FragLibraryTronDIATests::getMS2IonsTest() {

//    QSKIP("TEMP");

    ERR_INIT

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
    );
    QCOMPARE(e, eNoError);

    const PeptideSequenceChargeKey &peptideSequenceChargeKey = QStringLiteral("AGEVXVTAVAEHEK|3");

    QVector<MS2Ion> ms2Ions;
    e = fragLibraryTronDia.getMS2Ions(
            peptideSequenceChargeKey,
            &ms2Ions
            );

    QCOMPARE(peptideSequenceChargeKey, "AGEVXVTAVAEHEK|3");
    QCOMPARE(ms2Ions.size(), 43);

    const MS2Ion &ms2IonFront = ms2Ions.front();
    const MS2Ion &ms2IonBack = ms2Ions.back();

    qDebug() << ms2Ions;
    QCOMPARE(int(ms2IonFront.mz), 101);
    QCOMPARE(int(ms2IonFront.intensity * 10), 1);
    QCOMPARE(ms2IonFront.ionLabel, "a2");

    QCOMPARE(int(ms2IonBack.mz), int(1096.071));
    QCOMPARE(int(ms2IonBack.intensity * 10), 1);
    QCOMPARE(ms2IonBack.ionLabel, "y10");

}

void FragLibraryTronDIATests::getMS2IonsTest2() {

//    QSKIP("TEMP");

    ERR_INIT

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
    );
    QCOMPARE(e, eNoError);

    const PeptideSequenceChargeKey &peptideSequenceChargeKey = QStringLiteral("AGEVXVTAVAEHEK|3");

    const QPair<double, double> mzStartMzEnd = {101.072, 1096.599};

    QVector<MS2Ion> ms2Ions;
    e = fragLibraryTronDia.getMS2Ions(
            peptideSequenceChargeKey,
            mzStartMzEnd,
            &ms2Ions
    );

    QCOMPARE(peptideSequenceChargeKey, "AGEVXVTAVAEHEK|3");
    QCOMPARE(ms2Ions.size(), 41);

    const MS2Ion &ms2IonFront = ms2Ions.front();
    const MS2Ion &ms2IonBack = ms2Ions.back();

    qDebug() << ms2Ions;
    QCOMPARE(int(ms2IonFront.mz), 129);
    QCOMPARE(int(ms2IonFront.intensity * 100), 44);
    QCOMPARE(ms2IonFront.ionLabel, "b2");

    QCOMPARE(int(ms2IonBack.mz), 983);
    QCOMPARE(int(ms2IonBack.intensity * 100), 54);
    QCOMPARE(ms2IonBack.ionLabel, "y9");

}

void FragLibraryTronDIATests::getMS2IonsTest3() {

//    QSKIP("TEMP");

    ERR_INIT

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
    );
    QCOMPARE(e, eNoError);

    const PeptideSequenceChargeKey &peptideSequenceChargeKey = QStringLiteral("AGEVXVTAVAEHEK|3");

    const QPair<double, double> mzStartMzEnd = {101.072, 1096.599};

    const int topNIntense = 3;

    QVector<MS2Ion> ms2Ions;
    e = fragLibraryTronDia.getMS2Ions(
            peptideSequenceChargeKey,
            mzStartMzEnd,
            topNIntense,
            &ms2Ions
    );

    QCOMPARE(peptideSequenceChargeKey, "AGEVXVTAVAEHEK|3");
    QCOMPARE(ms2Ions.size(), topNIntense);

    const MS2Ion &ms2IonFront = ms2Ions.front();
    const MS2Ion &ms2IonBack = ms2Ions.back();

    qDebug() << ms2Ions;
    QCOMPARE(int(ms2IonFront.mz), 258);
    QCOMPARE(int(ms2IonFront.intensity * 100), 79);
    QCOMPARE(ms2IonFront.ionLabel, "b3");

    QCOMPARE(int(ms2IonBack.mz), 884);
    QCOMPARE(int(ms2IonBack.intensity), 1);
    QCOMPARE(ms2IonBack.ionLabel, "y8");

}

void FragLibraryTronDIATests::getMS2IonsTest4() {

    ERR_INIT

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            pythiaParameters(),
            fragFilePath()
    );
    QCOMPARE(e, eNoError);

    const PeptideSequenceChargeKey &peptideSequenceChargeKey = QStringLiteral("AGEVXVTAVAEHEK|3");

    const double mzTargetStart = 500;
    const double mzTargetEnd = 510;
    const int topNIntense = 3;

    QMap<PeptideStringWithMods, QVector<MS2Ion>> peptideStringWithModsVsMS2Ions;
    e = fragLibraryTronDia.getMS2Ions(
            mzTargetStart,
            mzTargetEnd,
            topNIntense,
            &peptideStringWithModsVsMS2Ions
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideStringWithModsVsMS2Ions.size(), 1960);

}


QTEST_MAIN(FragLibraryTronDIATests)

#include "FragLibraryTronDIATests.moc"
