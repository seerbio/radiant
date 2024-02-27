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
            100.0,
            4000.0,
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

    QString testFile = QStringLiteral("/home/anichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv");
//    testFile =

    const double massStart = 100.0;
    const double massEnd = 4000.0;

    FragLibTsvReader fragLibTsvReader;
    QVector<FragLibReaderRow> fragLibReaderTSVRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            massStart,
            massEnd,
            &fragLibReaderTSVRows
    );
    std::sort(fragLibReaderTSVRows.begin(), fragLibReaderTSVRows.end(), [](const FragLibReaderRow &l, const FragLibReaderRow &r){
        return l.peptideSequenceChargeKey < r.peptideSequenceChargeKey;
    });

    testFile = "/home/anichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF";
    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            testFile,
            massStart,
            massEnd,
            &fragLibReaderRows
            );
    std::sort(fragLibReaderRows.begin(), fragLibReaderRows.end(), [](const FragLibReaderRow &l, const FragLibReaderRow &r){
        return l.peptideSequenceChargeKey < r.peptideSequenceChargeKey;
    });

//    for (const FragLibReaderRow &flrr : fragLibReaderRows) {
//        qDebug() << flrr.peptideSequenceChargeKey << flrr.mzVals;
//    }

    QCOMPARE(fragLibReaderTSVRows.size(), fragLibReaderRows.size());

    for (int i = 0; i < fragLibReaderRows.size(); i++) {

        const FragLibReaderRow &tsv = fragLibReaderTSVRows.at(i);
        const FragLibReaderRow &reg = fragLibReaderRows.at(i);

        if (tsv.ionLabels != reg.ionLabels) {
            qDebug() << i;
            qDebug() << tsv.peptideSequenceChargeKey << reg.peptideSequenceChargeKey;
            qDebug() << tsv.mzVals << tsv.ionLabels;
            qDebug() << reg.mzVals << reg.ionLabels;
            qDebug() << "******";
        }

    }
}


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
