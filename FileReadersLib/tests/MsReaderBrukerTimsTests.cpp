//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MsReaderBrukerTims.h"


#include <QtTest/QtTest>

class MsReaderBrukerTimsTests : public QObject
{
    Q_OBJECT

public:
    MsReaderBrukerTimsTests() = default;
    ~MsReaderBrukerTimsTests() override = default;

private Q_SLOTS:
    void openFileTest();

};

void MsReaderBrukerTimsTests::openFileTest() {

    ERR_INIT

    const QString &openFileTest
        = QStringLiteral("/home/andrewnichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d");

    MsReaderBrukerTims reader;
    e = reader.openFile(openFileTest);
    QCOMPARE(e, eNoError);


}

QTEST_MAIN(MsReaderBrukerTimsTests)
#include "MsReaderBrukerTimsTests.moc"
