//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MsReaderBrukerTims.h"


#include <QtTest/QtTest>

#include "MsReaderMzMLMapped.h"

class MsReaderBrukerTimsTests : public QObject
{
    Q_OBJECT

public:
    MsReaderBrukerTimsTests() = default;
    ~MsReaderBrukerTimsTests() override = default;

private Q_SLOTS:
    void openFileTest();
    static void writeFrameTest();

};

void MsReaderBrukerTimsTests::writeFrameTest() {

    ERR_INIT

    QSKIP("Need to update fix writeFrame");

    const QString openFileTest
        = QStringLiteral("/home/andrewnichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d");

    constexpr float scanTime = 16.78975;

    MsReaderBrukerTims reader;
    e = reader.openFile(openFileTest);
    QCOMPARE(e, eNoError);

    e = reader.writeFrame(openFileTest + QString::number(scanTime) + ".prq", scanTime, 1);
    QCOMPARE(e, eNoError);
}

void MsReaderBrukerTimsTests::openFileTest() {

    // TODO make test

    ERR_INIT

    const QString &openFileTest
        = QStringLiteral("/home/andrewnichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d");

    MsReaderBrukerTims reader;
    e = reader.openFile(openFileTest);
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrame;
    e = reader.collateMS2MzTargetFrames(&diaTargetFrame);
    QCOMPARE(e, eNoError);

    MzTargetKey mzTargetKey = "409819";
    const QMap<ScanNumber, ScanPoints*> scanPoints = diaTargetFrame.value(mzTargetKey);

    QVector<MsScanInfo*> msScanInfos;
    e = reader.getMsScanInfos(mzTargetKey, &msScanInfos);
    QCOMPARE(e, eNoError);

    for (const MsScanInfo *msScanInfo : msScanInfos) {
        std::cout << msScanInfo->scanNumber << " " << msScanInfo->ionMobilityDriftTime << " " << msScanInfo->ionMobilityIndex << std::endl;
    }


}

QTEST_MAIN(MsReaderBrukerTimsTests)
#include "MsReaderBrukerTimsTests.moc"
