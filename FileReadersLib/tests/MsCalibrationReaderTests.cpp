//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsCalibrationReader.h"
#include "ParquetReader.h"

#include <QtConcurrent/QtConcurrent>
#include <QtTest/QtTest>

class MsCalibrationReaderTests : public QObject
{
    Q_OBJECT

public:
    MsCalibrationReaderTests() = default;
    ~MsCalibrationReaderTests() override = default;

private Q_SLOTS:

    void readWriteTest();

};

void MsCalibrationReaderTests::readWriteTest() {

    ERR_INIT

    MsCalibarationReaderRow row;
    row.peptideStringWithMods = QStringLiteral("KalliopeAndBellatrix");
    row.iRTPredicted = 66.6;
    row.scanTime = 66.7;
    row.scanNumber = 667;
    row.mzSearchedVec = {16.6, 16.7};
    row.mzFoundMeanVec = {26.6, 26.7};
    row.mzFoundStDevVec = {36.6, 36.7};
    row.intensityFoundMaxVec = {46.6, 46.7};

    const QVector<MsCalibarationReaderRow> rows = {row, row};

    const QString &msCalFilePath = QDir(qApp->applicationDirPath()).filePath("test.cal");

    e = ParquetReader::write(rows, msCalFilePath);
    QCOMPARE(e, eNoError);

    QVector<MsCalibarationReaderRow> rowsRead;
    e = ParquetReader::read(msCalFilePath, &rowsRead);
    QCOMPARE(e, eNoError);

    QCOMPARE(rowsRead.size(), rows.size());

    const MsCalibarationReaderRow &readRow = rowsRead.front();

    QCOMPARE(readRow.peptideStringWithMods, row.peptideStringWithMods);
    QCOMPARE(readRow.iRTPredicted, row.iRTPredicted);
    QCOMPARE(readRow.scanTime, row.scanTime);
    QCOMPARE(readRow.scanNumber, row.scanNumber);
    QCOMPARE(readRow.mzSearchedVec, row.mzSearchedVec);
    QCOMPARE(readRow.mzFoundMeanVec, row.mzFoundMeanVec);
    QCOMPARE(readRow.mzFoundStDevVec, row.mzFoundStDevVec);
    QCOMPARE(readRow.intensityFoundMaxVec, row.intensityFoundMaxVec);

}

QTEST_MAIN(MsCalibrationReaderTests)
#include "MsCalibrationReaderTests.moc"
