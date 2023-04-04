//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"


#include <QtTest/QtTest>

class TandemLibraryReaderTests : public QObject
{
    Q_OBJECT

public:
    TandemLibraryReaderTests() = default;
    ~TandemLibraryReaderTests() override = default;

private Q_SLOTS:

    void parseFastaFileTest();
};

void TandemLibraryReaderTests::parseFastaFileTest() {



}

QTEST_MAIN(TandemLibraryReaderTests)
#include "TandemLibraryReaderTests.moc"
