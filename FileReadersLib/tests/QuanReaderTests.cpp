//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "QuanReader.h"
#include "ParquetReader.h"

#include <QtConcurrent/QtConcurrent>
#include <QtTest/QtTest>

class QuanReaderTests : public QObject
{
    Q_OBJECT

public:
    QuanReaderTests() = default;
    ~QuanReaderTests() override = default;

private Q_SLOTS:

    static void readWrteTestCombined();

};

void QuanReaderTests::readWrteTestCombined() {

    QuanReaderRow qr;

    qr.peptideStringWithMods = "KALLIOPEANDBELLATRIX";
    qr.charge = 6;
    qr.targetKey = "666";
    qr.mzSearched1 = 666.1;
    qr.mzSearched2 = 666.2;
    qr.mzSearched3 = 666.3;
    qr.mzSearched4 = 666.4;
    qr.mzSearched5 = 666.5;
    qr.mzSearched6 = 666.6;
    qr.mzSearched7 = 666.7;
    qr.mzSearched8 = 666.8;
    qr.mzSearched9 = 666.9;
    qr.mzSearched10 = 666.10;
    qr.mzSearched11 = 666.11;
    qr.mzSearched12 = 666.12;
    qr.classifierScore = 0.00006;
    qr.discScore = 0.06;
    qr.qValue = 0.006;

    qr.intensityValsMz1 = {666.1, 666.2, 666.3, 666.4, 666.5, 666.6, 666.7};
    qr.intensityValsMz2 = {666.1, 666.2, 666.3, 666.4, 666.5, 666.6};
    qr.intensityValsMz3 = {666.1, 666.2, 666.3, 666.4, 666.5};
    qr.intensityValsMz4= {666.1, 666.2, 666.3, 666.4};
    qr.intensityValsMz5 = {666.1, 666.2, 666.3};
    qr.intensityValsMz6 = {666.1, 666.2};
    qr.intensityValsMz7 = {6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7};
    qr.intensityValsMz8  = {6.1, 6.2, 6.3, 6.4, 6.5, 6.6};
    qr.intensityValsMz9  = {6.1, 6.2, 6.3, 6.4, 6.5};
    qr.intensityValsMz10  = {6.1, 6.2, 6.3, 6.4};
    qr.intensityValsMz11  = {6.1, 6.2, 6.3};
    qr.intensityValsMz12 = {6.1, 6.2};

    qr.mzInterferences = {666.1};
    qr.isDecoy = true;

    ERR_INIT

    const QString testFilePath = QStringLiteral("test.txt");

    const QVector<QuanReaderRow> testData = {qr};
    e = ParquetReader::write(testData, testFilePath);
    QCOMPARE(e, eNoError);

    e = ErrorUtils::fileExists(testFilePath);
    QCOMPARE(e, eNoError);

    QVector<QuanReaderRow> readRows;
    e = ParquetReader::read(testFilePath, &readRows);
    QCOMPARE(e, eNoError);
    QCOMPARE(readRows.size(), 1);

    const QuanReaderRow readRow = readRows.at(0);

    QCOMPARE(readRow.peptideStringWithMods , qr.peptideStringWithMods);
    QCOMPARE(readRow.charge , qr.charge);
    QCOMPARE(readRow.targetKey , qr.targetKey);
    QCOMPARE(readRow.mzSearched1 , qr.mzSearched1);
    QCOMPARE(readRow.mzSearched2 , qr.mzSearched2);
    QCOMPARE(readRow.mzSearched3 , qr.mzSearched3);
    QCOMPARE(readRow.mzSearched4 , qr.mzSearched4);
    QCOMPARE(readRow.mzSearched5 , qr.mzSearched5);
    QCOMPARE(readRow.mzSearched6 , qr.mzSearched6);
    QCOMPARE(readRow.mzSearched7 , qr.mzSearched7);
    QCOMPARE(readRow.mzSearched8 , qr.mzSearched8);
    QCOMPARE(readRow.mzSearched9 , qr.mzSearched9);
    QCOMPARE(readRow.mzSearched10 , qr.mzSearched10);
    QCOMPARE(readRow.mzSearched11 , qr.mzSearched11);
    QCOMPARE(readRow.mzSearched12 , qr.mzSearched12);
    QCOMPARE(readRow.classifierScore , qr.classifierScore);
    QCOMPARE(readRow.discScore , qr.discScore);
    QCOMPARE(readRow.qValue , qr.qValue);
    QCOMPARE(readRow.intensityValsMz1 , qr.intensityValsMz1);
    QCOMPARE(readRow.intensityValsMz2 , qr.intensityValsMz2);
    QCOMPARE(readRow.intensityValsMz3 , qr.intensityValsMz3);
    QCOMPARE(readRow.intensityValsMz4 , qr.intensityValsMz4);
    QCOMPARE(readRow.intensityValsMz5 , qr.intensityValsMz5);
    QCOMPARE(readRow.intensityValsMz6 , qr.intensityValsMz6);
    QCOMPARE(readRow.intensityValsMz7 , qr.intensityValsMz7);
    QCOMPARE(readRow.intensityValsMz8 , qr.intensityValsMz8);
    QCOMPARE(readRow.intensityValsMz9 , qr.intensityValsMz9);
    QCOMPARE(readRow.intensityValsMz10 , qr.intensityValsMz10);
    QCOMPARE(readRow.intensityValsMz11 , qr.intensityValsMz11);
    QCOMPARE(readRow.intensityValsMz12 , qr.intensityValsMz12);
    QCOMPARE(readRow.mzInterferences , qr.mzInterferences);

    QCOMPARE(readRow.isDecoy , qr.isDecoy);

}

QTEST_MAIN(QuanReaderTests)
#include "QuanReaderTests.moc"
