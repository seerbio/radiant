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
    tpr.peptideSequenceChargeKey = "CHAUNCYANDFLOPS";
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

Err logic(const QString &testFilePath) {

    ERR_INIT

    const double massStart = 1000.0;
    const double massEnd = 1002.0;

    FragLibReader fragLibReader;
    e = fragLibReader.init(testFilePath); ree;


    QMap<PeptideStringWithMods, QVector<MS2Ion>> peptideStringWithModsVsMS2Ions;
    QMap<PeptideSequenceChargeKey, bool> peptideSequenceChargeKeyVsIsDecoy;
    QMap<PeptideSequenceChargeKey, double> peptideSequenceChargeKeyVsMass;
    e = fragLibReader.getMS2Ions(
            massStart,
            massEnd,
            &peptideStringWithModsVsMS2Ions,
            &peptideSequenceChargeKeyVsIsDecoy,
            &peptideSequenceChargeKeyVsMass
    ); ree;

    ERR_RETURN
}

void FragLibReaderTests::getSM2IonsTest() {

    const QString testFilePath
        = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";

    QVector<QString> paths(60, testFilePath);

    QFuture<Err> futures = QtConcurrent::mapped(
            paths,
            logic
            ) ;
    futures.waitForFinished();

}


QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
