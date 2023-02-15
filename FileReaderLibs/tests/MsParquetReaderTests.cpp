//
// Created by anichols on 11/07/2021.
//

#include "MsParquetReader.h"

#include <QtTest/QtTest>

class MsParquetReaderTests : public QObject
{
    Q_OBJECT

public:
    MsParquetReaderTests() = default;
    ~MsParquetReaderTests() override = default;

private Q_SLOTS:

    void checkParquetStatusTest();

};

void MsParquetReaderTests::checkParquetStatusTest() {

    const bool isOk = MsParquetReader::checkParquetStatus();
    QCOMPARE(isOk, true);

}


QTEST_MAIN(MsParquetReaderTests)
#include "MsParquetReaderTests.moc"
