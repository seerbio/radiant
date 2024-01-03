//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"

#include <QtConcurrent/QtConcurrent>
#include <QtTest/QtTest>

class FragLibReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibReaderTests() = default;
    ~FragLibReaderTests() override = default;

private Q_SLOTS:

    void writeTandemPredictionsAndReadTandemPredictionsCombinedTest();
    void getSM2IonsTest();

};

void FragLibReaderTests::writeTandemPredictionsAndReadTandemPredictionsCombinedTest() {

    FragLibReaderRow tpr;
    tpr.peptideSequenceChargeKey = "CHAUNCYANDFLOPS|666";
    tpr.intensityVals = {666.6, 66.6, 6.6};

    const QVector<FragLibReaderRow> tprs(10, tpr);

    ERR_INIT

    const QString &outputFilePath
            = QDir(qApp->applicationDirPath()).filePath("pred.tPreds");

    e = ParquetReader::write(
            tprs,
            outputFilePath
            );

    QCOMPARE(e, eNoError);

    QVector<FragLibReaderRow> readRows;
    e = ParquetReader::read(
            outputFilePath,
            &readRows
            );

    const FragLibReaderRow &tlrr = readRows.first();

    QCOMPARE(tlrr.peptideSequenceChargeKey, tpr.peptideSequenceChargeKey);
    QCOMPARE(tlrr.intensityVals, tpr.intensityVals);

}

void FragLibReaderTests::getSM2IonsTest() {

    ERR_INIT

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    QVector<FragLibReaderRow> fragLibReaderRows;

    e = FragLibReader::getFragLibReaderRows(
            testFilePath,
            100.0,
            40000.0,
            &fragLibReaderRows
            );

    QCOMPARE(e, eNoError);
    QCOMPARE(fragLibReaderRows.size(), 9581);

    e = FragLibReader::getFragLibReaderRows(
            testFilePath,
            900.0,
            600.0,
            &fragLibReaderRows
    );
    QCOMPARE(e, eError);

    e = FragLibReader::getFragLibReaderRows(
            testFilePath + "CHAUNCYANDFLOPS",
            600.0,
            900.0,
            &fragLibReaderRows
    );
    QCOMPARE(e, eFileError);

}

QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
