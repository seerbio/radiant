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
};

void FragLibTsvReaderTests::getFragLibReaderRowsTest() {

    ERR_INIT

    const QString &testFile = QDir(qApp->applicationDirPath()).filePath("tsv_library_trunc.tsv");

    FragLibTsvReader fragLibTsvReader;

    QVector<FragLibReaderRow> fragLibReaderRows;
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
    QVector<FragLibReaderRow> fragLibReaderTSVRows;
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


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
