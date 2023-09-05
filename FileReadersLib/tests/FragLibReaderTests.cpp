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
    void convertDIANNLibToFragLibTest();
};

void FragLibReaderTests::writeTandemPredictionsAndReadTandemPredictionsCombinedTest() {

    QSKIP("undo me baby");

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

    QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptide;
    e = FragLibReader::getMS2Ions(
            testFilePath,
            massStart,
            massEnd,
            100.0,
            2000.0,
            1000,
            &peptideSequenceChargeKeyVsCandidatePeptide
    ); ree;

    ERR_RETURN
}

void FragLibReaderTests::getSM2IonsTest() {

    QSKIP("undo me baby");

    const QString testFilePath
        = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";

    QVector<QString> paths(60, testFilePath);

    QFuture<Err> futures = QtConcurrent::mapped(
            paths,
            logic
            ) ;
    futures.waitForFinished();

}

void FragLibReaderTests::convertDIANNLibToFragLibTest() {

    ERR_INIT

    const QString diannLibraryFile = "/home/anichols/Desktop/Testing/LatestStuff/predicted_lib.predicted.speclib";

    e = FragLibReader::convertDIANNLibToFragLib(diannLibraryFile);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(FragLibReaderTests)
#include "FragLibReaderTests.moc"
