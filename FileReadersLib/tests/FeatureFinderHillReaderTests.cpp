//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FeatureFinderHillReader.h"
#include "ParquetReader.h"

#include <QtTest/QtTest>

class FeatureFinderHillReaderTests : public QObject
{
    Q_OBJECT

public:
    FeatureFinderHillReaderTests() = default;
    ~FeatureFinderHillReaderTests() override = default;

private Q_SLOTS:

    void readWriteTestCombined();
};

void FeatureFinderHillReaderTests::readWriteTestCombined() {

    ERR_INIT;

    const QString &testFilePath
        = QDir(qApp->applicationDirPath()).filePath("test.pythiaHills");

    const QVector<double> mzVals = {666.6, 66.6};
    const QVector<double> intensityVals = {0.666, 0.666};
    const QVector<int> scanNumbers = {666, 667};
    const QVector<int> scanNumberIndexes = {6666, 666};

    FeatureFinderHillReaderRow row;
    row.mzVals = mzVals;
    row.intensityVals = intensityVals;
    row.scanNumbers = scanNumbers;
    row.scanNumberIndexes = scanNumberIndexes;

    const QVector<FeatureFinderHillReaderRow> writeRows(6, row);

    e = ParquetReader::write(writeRows, testFilePath);
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHillReaderRow> readRows;
    e = ParquetReader::read(testFilePath, &readRows);
    QCOMPARE(e, eNoError);

    QCOMPARE(readRows.size(), writeRows.size());

    const FeatureFinderHillReaderRow &readRow = readRows.back();
    QCOMPARE(readRow.mzVals, mzVals);
    QCOMPARE(readRow.intensityVals, intensityVals);
    QCOMPARE(readRow.scanNumbers, scanNumbers);
    QCOMPARE(readRow.scanNumberIndexes, scanNumberIndexes);

}

QTEST_MAIN(FeatureFinderHillReaderTests)
#include "FeatureFinderHillReaderTests.moc"
