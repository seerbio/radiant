//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "CSVReader.h"

#include <QtTest/QtTest>

class CSVReaderTests : public QObject
{
    Q_OBJECT

public:
    CSVReaderTests() = default;
    ~CSVReaderTests() override = default;

private Q_SLOTS:
    void readTest();
    void writeTest();
};

void CSVReaderTests::readTest() {
    //TODO Implement
}

void CSVReaderTests::writeTest() {
    //TODO Implement
}

QTEST_MAIN(CSVReaderTests)
#include "CSVReaderTests.moc"
