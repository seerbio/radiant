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
};

void FragLibTsvReaderTests::getFragLibReaderRowsTest() {

    ERR_INIT

    //TODO make toy example and use proper pathing.
    const QString testFile = "/home/anichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv";

    FragLibTsvReader fragLibTsvReader;

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            100.0,
            4000.0,
            &fragLibReaderRows
            );

    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
