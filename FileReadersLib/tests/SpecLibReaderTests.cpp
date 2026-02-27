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
    static void peptideSearchTroubleShoot();

};


void SpecLibReaderTests::peptideSearchTroubleShoot() {

    QSKIP("Turned on only during troubleshooting");

    ERR_INIT

    const QString filePath = "/home/andrewnichols/Downloads/uniprot_homosapiens_isoforms_con_20221209.predicted(1).speclib";

    QList<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        filePath,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    const QString peptideSequence = "EAQEFVK";
    const int charge = 1;

    for (const FragLibReaderRow &flrr : fragLibReaderRows) {

        if (flrr.peptideSequenceChargeKey != peptideSequence + "|" + QString::number(charge)) {
            continue;
        }

        const QStringList ionLabels = flrr.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR);
        for (int i = 0; i < flrr.mzVals.size(); i++) {
            qDebug() << flrr.mzVals.at(i) << flrr.intensityVals.at(i) << ionLabels.at(i);
        }
        qDebug() << "Precursor charge" << charge;
    }

}

void SpecLibReaderTests::getFragLibReaerRowsTest() {

    QSKIP("Build a test for this when a small lib file is found");

    ERR_INIT

    // const QString specLibFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/lib.predicted.speclib");
    const QString specLibFile
        = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib");
    QList<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        specLibFile,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    const QString fragLibFFFile
    = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF");
    QList<FragLibReaderRow> fragLibReaderRowsFF;
    e = FragLibReader::getFragLibReaderRows(
        fragLibFFFile,
        false,
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
