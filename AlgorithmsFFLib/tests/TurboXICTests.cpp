//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "MsReaderMzML.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtTest/QtTest>

#include <algorithm>

class TurboXICTests : public QObject
{
    Q_OBJECT

public:
    TurboXICTests() = default;
    ~TurboXICTests() override = default;

private Q_SLOTS:

    void initTest();
    void extractPointsTest();
    void scanRestrictedQueryPreservesOrderAndValues();
    void turboXICUtility();


private:

    QMap<int, QVector<PointFF>> buildPoints() const;

};


QMap<int, QVector<PointFF>> TurboXICTests::buildPoints() const {

    QMap<int, QVector<PointFF>> points;
    points.insert(1, {PointFF(100.11, 100.1), PointFF(200.11, 200.1)});
    points.insert(2, {PointFF(100.12, 100.2), PointFF(200.12, 200.2)});
    points.insert(3, {PointFF(100.13, 100.3), PointFF(200.13, 200.3)});
    points.insert(4, {PointFF(100.14, 100.4), PointFF(200.14, 200.4)});
    points.insert(5, {PointFF(100.15, 100.5), PointFF(100.151, 100.5), PointFF(200.15, 200.5)});


    return points;
}

void TurboXICTests::initTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanPoints*> emptyPoints;
    e = turboXIC.init(emptyPoints);
    QCOMPARE(e, eEmptyContainerError);

}

void TurboXICTests::extractPointsTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    //TODO fix test

    const XICPoints xicPoints = turboXIC.extractPointsXIC(100.0, 100.13);

    QCOMPARE(xicPoints.size(), 3);
    QCOMPARE(xicPoints.front().scanNumber, 1);
    QVERIFY(MathUtils::tSame(xicPoints.front().intensity, 100.1f));
    QVERIFY(MathUtils::tSame(xicPoints.back().intensity, 100.3f));

}

void TurboXICTests::scanRestrictedQueryPreservesOrderAndValues() {
    QMap<int, QVector<PointFF>> storage = buildPoints();
    QMap<ScanNumber, ScanPoints*> points;
    for (auto it = storage.begin(); it != storage.end(); ++it)
        points.insert(it.key(), &it.value());

    TurboXIC turboXIC;
    QCOMPARE(turboXIC.init(points), eNoError);
    const auto all = turboXIC.extractPointsXIC(100.0f, 100.2f);
    const auto compareRestricted = [&](ScanNumber minScan, ScanNumber maxScan,
                                       bool includeEndpoints) {
        auto expected = all;
        expected.erase(std::remove_if(expected.begin(), expected.end(),
            [&](const XICPoint &point) {
                return includeEndpoints
                    ? !(minScan <= point.scanNumber && point.scanNumber <= maxScan)
                    : !(minScan < point.scanNumber && point.scanNumber < maxScan);
            }), expected.end());
        const auto actual = turboXIC.extractPointsXIC(
            100.0f, 100.2f, minScan, maxScan, includeEndpoints);
        QCOMPARE(actual.size(), expected.size());
        for (size_t i = 0; i < actual.size(); ++i) {
            QCOMPARE(actual[i].scanNumber, expected[i].scanNumber);
            QCOMPARE(actual[i].mz, expected[i].mz);
            QCOMPARE(actual[i].intensity, expected[i].intensity);
        }
    };
    compareRestricted(1, 3, true);
    compareRestricted(1, 3, false);
    compareRestricted(4, 4, true);
    compareRestricted(4, 4, false);
}

void TurboXICTests::turboXICUtility() {

    ERR_INIT

    QSKIP("uncomment for troubleshooting");

//    const QString msDataFilePath
//            = QStringLiteral("/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");
//
//    const MzTargetKey uniqueMsInfoScanKey = "454957";
//    double target = 454.957;
//    double window = 5.5;
//
//    MsFrame msFrame;
//    e = MsFrame::buildMsFrame(
//            msDataFilePath,
//            uniqueMsInfoScanKey,
//            {target - window, target + window},
//            &msFrame
//    );
//    QCOMPARE(e, eNoError);
//
//    e = msFrame.smoothFrame(
//            3,
//            1.0,
//            2,
//            1500.0
//            );
//    QCOMPARE(e, eNoError);
//
//    const QMap<int, QVector<PointFF>> &points = msFrame.scanNumberVsScanPoints();
//
//    TurboXIC turboXIC;
//    e = turboXIC.init(points);
//    QCOMPARE(e, eNoError);
//
//    const double mzCenter = 523.26232347;
//    const double ppmTol = 50;
//    const double massTol = MathUtils::calculatePPM(mzCenter, ppmTol);
//
//    const XICPoints xicPoints = turboXIC.extractPointsXIC(
//            mzCenter - massTol,
//            mzCenter + massTol,
//            0,
//            26000
//            );
//
//    qDebug() << xicPoints.scanNumbersVsIntensityVals.size();
//    const QVector<PointFF> vec = ParallelUtils::convertMapToPoints(xicPoints.scanNumbersVsIntensityVals);
//
//    for (const PointFF &p : vec) {
//        qDebug() << p;
//    }
//
//    e = MsUtils::writePointsToCSV(vec, "xic.csv");
//    QCOMPARE(e, eNoError);

}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
