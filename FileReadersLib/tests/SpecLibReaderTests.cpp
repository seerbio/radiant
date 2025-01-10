//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReaderRow.h"
#include "SpecLibReader.h"

#include <QtTest/QtTest>

#include "FragLibReader.h"

class SpecLibReaderTests : public QObject
{
    Q_OBJECT

public:
    SpecLibReaderTests() = default;
    ~SpecLibReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaerRowsTest();

};


void SpecLibReaderTests::getFragLibReaerRowsTest() {

    // QSKIP("Build a test for this when a small lib file is found");

    ERR_INIT

    // const QString specLibFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/lib.predicted.speclib");
    const QString specLibFile
        = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib");
    QVector<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        specLibFile,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    const QString fragLibFFFile
    = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF");
    QVector<FragLibReaderRow> fragLibReaderRowsFF;
    e = FragLibReader::getFragLibReaderRows(
        fragLibFFFile,
        &fragLibReaderRowsFF
        );
    QCOMPARE(e, eNoError);

    qDebug() << fragLibReaderRowsFF.size() << fragLibReaderRows.size() << "Sizes";

    for (int i = 0; i < fragLibReaderRows.size(); i++) {
        QCOMPARE(fragLibReaderRows.at(i).peptideSequenceChargeKey, fragLibReaderRowsFF.at(i).peptideSequenceChargeKey);
        QCOMPARE(fragLibReaderRows.at(i).intensityVals, fragLibReaderRowsFF.at(i).intensityVals);
        // if (fragLibReaderRows.at(i).ionLabels != fragLibReaderRowsFF.at(i).ionLabels) {
        //     qDebug() << fragLibReaderRows.at(i).ionLabels << fragLibReaderRowsFF.at(i).ionLabels;
        //     qDebug() << fragLibReaderRows.at(i).intensityVals << fragLibReaderRowsFF.at(i).intensityVals;
        //
        // }


    }


}


QTEST_MAIN(SpecLibReaderTests)
#include "SpecLibReaderTests.moc"
