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

QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
