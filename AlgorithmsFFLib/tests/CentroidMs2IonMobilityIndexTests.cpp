//
// Created by Codex on 7/29/26.
//

#include "CentroidMs2IonMobilityIndex.h"

#include "ErrorUtils.h"
#include "MsFrame.h"

#include <QtTest/QtTest>

class CentroidMs2IonMobilityIndexTests : public QObject
{
    Q_OBJECT

public:
    CentroidMs2IonMobilityIndexTests() = default;
    ~CentroidMs2IonMobilityIndexTests() override = default;

private Q_SLOTS:
    static void extractFiltersMzFrameAndIonMobilityTest();
};

void CentroidMs2IonMobilityIndexTests::extractFiltersMzFrameAndIonMobilityTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    scanNumberVsScanPoints.insert(100, {{500.000f, 100.0f}, {501.000f, 999.0f}});
    scanNumberVsScanPoints.insert(200, {{500.001f, 200.0f}});
    scanNumberVsScanPoints.insert(300, {{500.002f, 300.0f}});

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

    QMap<ScanNumber, const TimsbukAlignedPointData*> scanNumberVsAlignedPointData;
    const TimsbukAlignedPointData alignedScan100{{1.00f, 1.15f}};
    const TimsbukAlignedPointData alignedScan200{{1.00f}};
    const TimsbukAlignedPointData alignedScan300{{1.20f}};
    scanNumberVsAlignedPointData.insert(100, &alignedScan100);
    scanNumberVsAlignedPointData.insert(200, &alignedScan200);
    scanNumberVsAlignedPointData.insert(300, &alignedScan300);

    CentroidMs2IonMobilityIndex index;
    e = index.init(
        scanNumberVsScanPoints,
        scanNumberVsAlignedPointData,
        msFrame
        );
    QCOMPARE(e, eNoError);
    QCOMPARE(index.pointCount(), 4);

    float driftTime = -1.0f;
    QVERIFY(index.driftTimeFromIonMobilityIndex(10000, &driftTime));
    QVERIFY(MathUtils::tSame(driftTime, 1.00f));
    QVERIFY(index.driftTimeFromIonMobilityIndex(12000, &driftTime));
    QVERIFY(MathUtils::tSame(driftTime, 1.20f));
    QVERIFY(!index.driftTimeFromIonMobilityIndex(13000, &driftTime));

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
    QCOMPARE(xicPoints.at(0).ionMobilityIndex, 10000);
    QCOMPARE(xicPoints.at(1).ionMobilityIndex, 10000);
    QVERIFY(MathUtils::tSame(xicPoints.at(0).intensity, 100.0f));
    QVERIFY(MathUtils::tSame(xicPoints.at(1).intensity, 200.0f));

    const XICPoints narrowMobilityPoints = index.extractPointsXIC(
        500.99f,
        501.01f,
        -1,
        2,
        1.10f,
        1.16f
        );

    QCOMPARE(narrowMobilityPoints.size(), 1);
    QCOMPARE(narrowMobilityPoints.at(0).scanNumber, 0);
    QCOMPARE(narrowMobilityPoints.at(0).ionMobilityIndex, 11500);
    QVERIFY(MathUtils::tSame(narrowMobilityPoints.at(0).intensity, 999.0f));
}

QTEST_MAIN(CentroidMs2IonMobilityIndexTests)
#include "CentroidMs2IonMobilityIndexTests.moc"
