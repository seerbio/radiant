//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "FragLibTsvReader.h"


#include <QFile>
#include <QTextStream>
#include <QtTest/QtTest>

class FragLibTsvReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibTsvReaderTests() = default;
    ~FragLibTsvReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaderRowsTest();
    static void getFragLibReaderRowsFiltersInvalidSequencesTest();

    static void compareTest();
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

void FragLibTsvReaderTests::getFragLibReaderRowsFiltersInvalidSequencesTest() {
    ERR_INIT

    const QString testFile = QDir(qApp->applicationDirPath()).filePath("FragLibTsvReaderTests.invalid_sequences.tsv");
    QFile file(testFile);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text));

    QTextStream out(&file);
    out << "PrecursorMz\tProductMz\tTr_recalibrated\tIonMobility\tLibraryIntensity\tdecoy\tModifiedPeptide\tPrecursorCharge\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\n";
    out << "500.0\t100.0\t10.0\t1.0\t1000.0\t0\tACDE\t2\tb\t1\t1\n";
    out << "500.0\t200.0\t10.0\t1.0\t500.0\t0\tACDE\t2\ty\t1\t1\n";
    out << "600.0\t110.0\t20.0\t1.1\t900.0\t0\tACDX\t2\tb\t1\t1\n";
    out << "600.0\t210.0\t20.0\t1.1\t450.0\t0\tACDX\t2\ty\t1\t1\n";
    file.close();

    QList<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(testFile, &fragLibReaderRows);
    QCOMPARE(e, eNoError);

    QCOMPARE(fragLibReaderRows.size(), 1);
    QCOMPARE(fragLibReaderRows.front().peptideSequenceChargeKey, QString("ACDE|2"));
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


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
