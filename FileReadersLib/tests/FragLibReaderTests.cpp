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

    void readWrteTestCombined();
    void getSM2IonsTest();
    void filterInvalidSequencesFragLibFFTest();

};

void FragLibReaderTests::readWrteTestCombined() {

    FragLibReaderRow tpr;
    tpr.peptideSequenceChargeKey = "CHAUNCYANDFLOPS|666";
    tpr.intensityVals = {666.6, 66.6, 6.6};

    QList<FragLibReaderRow> tprs;
    for (int i = 0; i < 10; i++) {
        tprs.push_back(tpr);
    }

    ERR_INIT

    const QString &outputFilePath
            = QDir(qApp->applicationDirPath()).filePath("pred.tPreds");

    e = ParquetReader::write(
            tprs,
            outputFilePath
            );

    QCOMPARE(e, eNoError);

    QList<FragLibReaderRow> readRows;
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

    QList<FragLibReaderRow> fragLibReaderRows;

    e = FragLibReader::getFragLibReaderRows(
            testFilePath,
            false,
            &fragLibReaderRows
            );

    QCOMPARE(e, eNoError);
    QCOMPARE(fragLibReaderRows.size(), 9581);

    e = FragLibReader::getFragLibReaderRows(
            testFilePath + "CHAUNCYANDFLOPS",
            false,
            &fragLibReaderRows
    );
    QCOMPARE(e, eFileError);

}

void FragLibReaderTests::filterInvalidSequencesFragLibFFTest() {
    ERR_INIT

    FragLibReaderRow valid;
    valid.peptideSequenceChargeKey = "ACDEFGHIK|2";
    valid.precursorCharge = 2;
    valid.iRT = 1.0f;
    valid.intensityVals = {1.0f, 0.5f};
    valid.mzVals = {100.0f, 200.0f};
    valid.ionLabels = "b1;y1";
    valid.mass = 1000.0f;

    FragLibReaderRow invalid = valid;
    invalid.peptideSequenceChargeKey = "ACDX|2";

    QList<FragLibReaderRow> rowsIn = {valid, invalid};

    const QString filePath = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.invalid_sequences.fragLibFF");
    e = ParquetReader::write(rowsIn, filePath);
    QCOMPARE(e, eNoError);

    QList<FragLibReaderRow> rowsOut;
    e = FragLibReader::getFragLibReaderRows(filePath, false, &rowsOut);
    QCOMPARE(e, eNoError);

    QCOMPARE(rowsOut.size(), 1);
    QCOMPARE(rowsOut.front().peptideSequenceChargeKey, valid.peptideSequenceChargeKey);
}

QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
