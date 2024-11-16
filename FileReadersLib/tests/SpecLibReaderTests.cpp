//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReaderRow.h"
#include "SpecLibReader.h"

#include <QtTest/QtTest>

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

    ERR_INIT

    const QString specLibFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/lib.predicted.speclib");

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        specLibFile,
        &fragLibReaderRows
        );

    QCOMPARE(e, eNoError);
}



QTEST_MAIN(SpecLibReaderTests)
#include "SpecLibReaderTests.moc"
