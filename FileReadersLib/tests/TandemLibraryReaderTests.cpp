//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "TandemLibraryReader.h"


#include <QtTest/QtTest>

class TandemLibraryReaderTests : public QObject
{
    Q_OBJECT

public:
    TandemLibraryReaderTests() = default;
    ~TandemLibraryReaderTests() override = default;

private Q_SLOTS:

    void writeTandemPredictionsAndReadTandemPredictionsCombinedTest();
};

void TandemLibraryReaderTests::writeTandemPredictionsAndReadTandemPredictionsCombinedTest() {

    TandemLibraryReaderRow tpr;
    tpr.peptideSequenceChargeKey = "CHAUNCYANDFLOPS";
    tpr.intensityVals = {666.6, 66.6, 6.6};
    tpr.ionLabels = QStringList({"a", "b", "c"});

    const QVector<TandemLibraryReaderRow> tprs(10, tpr);

    ERR_INIT

    const QString &outputFilePath
            = QDir(qApp->applicationDirPath()).filePath("pred.tPreds");

    e = TandemLibraryReader::writeTandemPredictions(
            tprs,
            outputFilePath
            );

    QCOMPARE(e, eNoError);

    QVector<TandemLibraryReaderRow> readRows;
    e = TandemLibraryReader::readTandemPredictions(
            outputFilePath.toStdString(),
            &readRows
            );

    const TandemLibraryReaderRow &tlrr = readRows.first();

    QCOMPARE(tlrr.peptideSequenceChargeKey, tpr.peptideSequenceChargeKey);
    QCOMPARE(tlrr.intensityVals, tpr.intensityVals);
    QCOMPARE(tlrr.ionLabels, tpr.ionLabels);
}

QTEST_MAIN(TandemLibraryReaderTests)
#include "TandemLibraryReaderTests.moc"
