//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtTest/QtTest>

class TurboXICTests : public QObject
{
    Q_OBJECT

public:
    TurboXICTests() = default;
    ~TurboXICTests() override = default;

private Q_SLOTS:

    void initTest();
    void extractPointsTest();

    void turboXICUtility();


private:

    QMap<int, QVector<PointDF>> buildPoints() const;

};


QMap<int, QVector<PointDF>> TurboXICTests::buildPoints() const {

    QMap<int, QVector<PointDF>> points;
    points.insert(1, {PointDF(100.11, 100.1), PointDF(200.11, 200.1)});
    points.insert(2, {PointDF(100.12, 100.2), PointDF(200.12, 200.2)});
    points.insert(3, {PointDF(100.13, 100.3), PointDF(200.13, 200.3)});
    points.insert(4, {PointDF(100.14, 100.4), PointDF(200.14, 200.4)});
    points.insert(5, {PointDF(100.15, 100.5), PointDF(100.151, 100.5), PointDF(200.15, 200.5)});


    return points;
}


void TurboXICTests::initTest() {

    QMap<int, QVector<PointDF>> _points = buildPoints();

    QMap<int, QVector<PointDF>*> points;
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

    QMap<int, QVector<PointDF>> _points = buildPoints();

    QMap<int, QVector<PointDF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    const XICPoints xicPoints = turboXIC.extractPointsXIC(100.0, 100.13, 2, 4);
    QCOMPARE(xicPoints.scanNumbersVsIntensity.size(), 2);
    QCOMPARE(xicPoints.scanNumbersVsIntensity.firstKey(), 2);
    QCOMPARE(xicPoints.scanNumbersVsIntensity.first(), 100.2);
    QCOMPARE(xicPoints.scanNumbersVsIntensity.lastKey(), 3);
    QCOMPARE(xicPoints.scanNumbersVsIntensity.last(), 100.3);

    const XICPoints xicPointsSum = turboXIC.extractPointsXIC(100.0, 100.17, 2, 5);
    QCOMPARE(xicPointsSum.scanNumbersVsIntensity.size(), 4);
    QCOMPARE(xicPointsSum.scanNumbersVsIntensity.last(), 201.);

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
//    const QMap<int, QVector<PointDF>> &points = msFrame.scanNumberVsScanPoints();
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
//    const QVector<PointDF> vec = ParallelUtils::convertMapToPoints(xicPoints.scanNumbersVsIntensityVals);
//
//    for (const PointDF &p : vec) {
//        qDebug() << p;
//    }
//
//    e = MsUtils::writePointsToCSV(vec, "xic.csv");
//    QCOMPARE(e, eNoError);

}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
