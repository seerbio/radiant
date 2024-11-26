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
    qr.scanTime = 6.6;
    qr.scanTimeStart = 6.5;
    qr.scanTimeEnd = 6.7;
    qr.intensityValMz1 = 666.1;
    qr.intensityValMz2 = 666.2;
    qr.intensityValMz3 = 666.3;
    qr.intensityValMz4= 666.4;
    qr.intensityValMz5 = 666.5;
    qr.intensityValMz6 = 666.6;

    qr.cosineSimSum100 = 666.777;

    qr.isDecoy = true;
    qr.scanTimePredicted = 8.6;


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
    QCOMPARE(readRow.scanTime , qr.scanTime);
    QCOMPARE(readRow.scanTimeStart , qr.scanTimeStart);
    QCOMPARE(readRow.scanTimeEnd , qr.scanTimeEnd);
    QCOMPARE(readRow.intensityValMz1 , qr.intensityValMz1);
    QCOMPARE(readRow.intensityValMz2 , qr.intensityValMz2);
    QCOMPARE(readRow.intensityValMz3 , qr.intensityValMz3);
    QCOMPARE(readRow.intensityValMz4 , qr.intensityValMz4);
    QCOMPARE(readRow.intensityValMz5 , qr.intensityValMz5);
    QCOMPARE(readRow.intensityValMz6 , qr.intensityValMz6);
    QCOMPARE(readRow.cosineSimSum100 , qr.cosineSimSum100);
    QCOMPARE(readRow.isDecoy , qr.isDecoy);
    QCOMPARE(readRow.scanTimePredicted , qr.scanTimePredicted);

}

QTEST_MAIN(QuanReaderTests)
#include "QuanReaderTests.moc"
