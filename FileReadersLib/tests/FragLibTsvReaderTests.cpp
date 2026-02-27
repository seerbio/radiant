//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "FragLibTsvReader.h"


#include <QtTest/QtTest>

class FragLibTsvReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibTsvReaderTests() = default;
    ~FragLibTsvReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaderRowsTest();

    static void compareTest();

    static void inferIonLabelsModeSwitchTest();
};

void FragLibTsvReaderTests::getFragLibReaderRowsTest() {

    ERR_INIT

    const QString &testFile = QDir(qApp->applicationDirPath()).filePath("tsv_library_trunc.tsv");

    FragLibTsvReader fragLibTsvReader;

    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(fragLibReaderRows.size(), 57);

    const FragLibReaderRow &fragLibReaderRow = fragLibReaderRows.at(47);
    qDebug() << fragLibReaderRow.peptideSequenceChargeKey
             << fragLibReaderRow.precursorCharge
             << fragLibReaderRow.mzVals
             << fragLibReaderRow.isDecoy
             << fragLibReaderRow.intensityVals
             << fragLibReaderRow.mass
             << fragLibReaderRow.iRT
             << fragLibReaderRow.ionLabels;

    QCOMPARE(fragLibReaderRow.precursorCharge, 1);
    QCOMPARE(fragLibReaderRow.mzVals.size(), 7);
    QCOMPARE(fragLibReaderRow.intensityVals.size(), 7);
    QCOMPARE(fragLibReaderRow.isDecoy, 0);
    QCOMPARE(static_cast<int>(fragLibReaderRow.mass), 600);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iRT), -38);
    QCOMPARE(fragLibReaderRow.ionLabels.size(), 20);

    QCOMPARE(e, eNoError);

}

void FragLibTsvReaderTests::compareTest() {

    QSKIP("For testing");

    ERR_INIT

    QString testFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv");
    // testFile = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib_trunc.tsv";
    // testFile = "/home/andrewnichols/Desktop/Data/Libraries/Test.tsv";
    // testFile = "/home/andrewnichols/Downloads/pythia-firstpass-lib.tsv";
    testFile = "/home/andrewnichols/Downloads/lib-mbr-2025-02-06.tsv";

    FragLibTsvReader fragLibTsvReader;
    QList<FragLibReaderRow> fragLibReaderTSVRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderTSVRows
    );


    // for (int i = 0; i < fragLibReaderTSVRows.size(); i++) {
    //
    //     const FragLibReaderRow &tsv = fragLibReaderTSVRows.at(i);
    //     qDebug()
    //     << i
    //     << tsv.peptideSequenceChargeKey
    //     << tsv.mzVals
    //     << tsv.intensityVals
    //     << tsv.ionLabels;
    //
    // }
}

void FragLibTsvReaderTests::inferIonLabelsModeSwitchTest() {

    const PeptideStringWithMods peptideStringWithMods("VTWLSPTNK");
    const bool isDecoy = true;
    const int charge = 1;

    QString ionLabels;

    // Default annotation mode (no decoy shift applied).
    {
        const QVector<float> mzVals = {
                546.28800f, // y5
                587.31900f  // b5
        };
        Err e = FragLibTsvReader::inferIonLabelsForTest(
                mzVals,
                peptideStringWithMods,
                isDecoy,
                charge,
                false,
                &ionLabels
                );
        QCOMPARE(e, eNoError);
        QCOMPARE(ionLabels, QString("y5;b5"));
    }

    // Terminal-by-penultimate mode: y7 and b7 retain terminal-series masses.
    {
        const QVector<float> mzVals = {
                560.30400f, // y5
                573.30300f  // b5
        };
        Err e = FragLibTsvReader::inferIonLabelsForTest(
                mzVals,
                peptideStringWithMods,
                isDecoy,
                charge,
                true,
                &ionLabels
                );
        QCOMPARE(e, eNoError);
        QCOMPARE(ionLabels, QString("y5;b5"));
    }
}


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
