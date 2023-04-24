//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"


#include <QtTest/QtTest>

class FragLibReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibReaderTests() = default;
    ~FragLibReaderTests() override = default;

private Q_SLOTS:

    void writeTandemPredictionsAndReadTandemPredictionsCombinedTest();
};

void FragLibReaderTests::writeTandemPredictionsAndReadTandemPredictionsCombinedTest() {

    FragLibReaderRow tpr;
    tpr.peptideSequenceChargeKey = "CHAUNCYANDFLOPS";
    tpr.intensityVals = {666.6, 66.6, 6.6};
//    tpr.ionLabels = QStringList({"a", "b", "c"});

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
//    QCOMPARE(tlrr.ionLabels, tpr.ionLabels);
}

QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
