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
    qr.classifierScore = 0.00006;
    qr.discScore = 0.06;
    qr.qValue = 0.006;

    qr.scanTimesMz1 = {6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7};
    qr.intensityValsMz1 = {666.1, 666.2, 666.3, 666.4, 666.5, 666.6, 666.7};
    qr.scanTimesMz2  = {6.1, 6.2, 6.3, 6.4, 6.5, 6.6};
    qr.intensityValsMz2 = {666.1, 666.2, 666.3, 666.4, 666.5, 666.6};
    qr.scanTimesMz3  = {6.1, 6.2, 6.3, 6.4, 6.5};
    qr.intensityValsMz3 = {666.1, 666.2, 666.3, 666.4, 666.5};
    qr.scanTimesMz4  = {6.1, 6.2, 6.3, 6.4};
    qr.intensityValsMz4= {666.1, 666.2, 666.3, 666.4};
    qr.scanTimesMz5  = {6.1, 6.2, 6.3};
    qr.intensityValsMz5 = {666.1, 666.2, 666.3};
    qr.scanTimesMz6 = {6.1, 6.2};
    qr.intensityValsMz6 = {666.1, 666.2};

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
    QCOMPARE(readRow.classifierScore , qr.classifierScore);
    QCOMPARE(readRow.discScore , qr.discScore);
    QCOMPARE(readRow.qValue , qr.qValue);
    QCOMPARE(readRow.scanTimesMz1 , qr.scanTimesMz1);
    QCOMPARE(readRow.intensityValsMz1 , qr.intensityValsMz1);
    QCOMPARE(readRow.scanTimesMz2 , qr.scanTimesMz2);
    QCOMPARE(readRow.intensityValsMz2 , qr.intensityValsMz2);
    QCOMPARE(readRow.scanTimesMz3 , qr.scanTimesMz3);
    QCOMPARE(readRow.intensityValsMz3 , qr.intensityValsMz3);
    QCOMPARE(readRow.scanTimesMz4 , qr.scanTimesMz4);
    QCOMPARE(readRow.intensityValsMz4 , qr.intensityValsMz4);
    QCOMPARE(readRow.scanTimesMz5 , qr.scanTimesMz5);
    QCOMPARE(readRow.intensityValsMz5 , qr.intensityValsMz5);
    QCOMPARE(readRow.scanTimesMz6 , qr.scanTimesMz6);
    QCOMPARE(readRow.intensityValsMz6 , qr.intensityValsMz6);

}

QTEST_MAIN(QuanReaderTests)
#include "QuanReaderTests.moc"
