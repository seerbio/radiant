//
// Created by Codex on 6/20/26.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "TimsMs2IonMobilityIndex.h"

#include <QtTest/QtTest>

class TimsMs2IonMobilityIndexTests : public QObject
{
    Q_OBJECT

public:
    TimsMs2IonMobilityIndexTests() = default;
    ~TimsMs2IonMobilityIndexTests() override = default;

private Q_SLOTS:
    static void extractFiltersMzFrameAndIonMobilityTest();
};

void TimsMs2IonMobilityIndexTests::extractFiltersMzFrameAndIonMobilityTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    scanNumberVsScanPoints.insert(100, {{500.0f, 1.0f}});
    scanNumberVsScanPoints.insert(200, {{500.0f, 1.0f}});
    scanNumberVsScanPoints.insert(300, {{500.0f, 1.0f}});

    QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs;
    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); ++it) {
        scanNumberVsScanPointsPntrs.insert(it.key(), &it.value());
    }

    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    scanNumberVsScanTime.insert(100, 10.0f);
    scanNumberVsScanTime.insert(200, 20.0f);
    scanNumberVsScanTime.insert(300, 30.0f);

    MsFrame msFrame;
    e = msFrame.init(scanNumberVsScanPointsPntrs, scanNumberVsScanTime);
    QCOMPARE(e, eNoError);

    QMap<FrameNumberTIMS, Ms2FrameTIMS> frameNumberVsMs2FrameTims;
    frameNumberVsMs2FrameTims[100][10] = {{500.000f, 100.0f}, {501.000f, 999.0f}};
    frameNumberVsMs2FrameTims[200][10] = {{500.001f, 200.0f}};
    frameNumberVsMs2FrameTims[200][20] = {{500.002f, 300.0f}};
    frameNumberVsMs2FrameTims[300][10] = {{500.003f, 400.0f}};

    QMap<FrameIndex, double> ionMobilityIndexVsDriftTime;
    ionMobilityIndexVsDriftTime.insert(10, 1.00);
    ionMobilityIndexVsDriftTime.insert(20, 1.20);
    ionMobilityIndexVsDriftTime.insert(30, 1.40);

    TimsMs2IonMobilityIndex index;
    e = index.init(
        frameNumberVsMs2FrameTims,
        msFrame,
        ionMobilityIndexVsDriftTime
        );
    QCOMPARE(e, eNoError);
    QCOMPARE(index.pointCount(), 5);

    float driftTime = -1.0f;
    QVERIFY(index.driftTimeFromIonMobilityIndex(20, &driftTime));
    QVERIFY(MathUtils::tSame(driftTime, 1.20f));
    QVERIFY(index.driftTimeFromIonMobilityIndex(30, &driftTime));
    QVERIFY(MathUtils::tSame(driftTime, 1.40f));
    QVERIFY(!index.driftTimeFromIonMobilityIndex(40, &driftTime));

    const XICPoints xicPoints = index.extractPointsXIC(
        499.99f,
        500.01f,
        -1,
        2,
        0.95f,
        1.05f
        );

    QCOMPARE(xicPoints.size(), 2);
    QCOMPARE(xicPoints.at(0).scanNumber, 0);
    QCOMPARE(xicPoints.at(1).scanNumber, 1);
    QCOMPARE(xicPoints.at(0).ionMobilityIndex, 10);
    QCOMPARE(xicPoints.at(1).ionMobilityIndex, 10);
    QVERIFY(MathUtils::tSame(xicPoints.at(0).intensity, 100.0f));
    QVERIFY(MathUtils::tSame(xicPoints.at(1).intensity, 200.0f));

    QMap<IonMobilityIndex, double> mobilityProfile;
    float apexIntensity = -1.0f;
    float apexDeltaAbs = -1.0f;
    QVERIFY(index.extractMobilityProfile(
        499.99f,
        500.01f,
        -1,
        2,
        0.95f,
        1.05f,
        1.02f,
        &mobilityProfile,
        &apexIntensity,
        &apexDeltaAbs
        ));

    QCOMPARE(mobilityProfile.size(), 1);
    QVERIFY(MathUtils::tSame(static_cast<float>(mobilityProfile.value(10)), 300.0f));
    QVERIFY(MathUtils::tSame(apexIntensity, 200.0f));
    QVERIFY(MathUtils::tSame(apexDeltaAbs, 0.02f));
}

QTEST_MAIN(TimsMs2IonMobilityIndexTests)
#include "TimsMs2IonMobilityIndexTests.moc"
