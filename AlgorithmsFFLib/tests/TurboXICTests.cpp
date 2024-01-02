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
    void mzMLMemoryTest();


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

//    const XICPoints xicPoints = turboXIC.extractPointsXIC(100.0, 100.13, 2, 4);
//    QCOMPARE(xicPoints.size(), 2);
//    QCOMPARE(xicPoints.front().scanNumber, 2);
//    QCOMPARE(MathUtils::pRound(xicPoints.scanNumbersVsIntensity.first(), 1), 100.2);
//    QCOMPARE(xicPoints.scanNumbersVsIntensity.lastKey(), 3);
//    QCOMPARE(MathUtils::pRound(xicPoints.scanNumbersVsIntensity.last(), 1), 100.3);
//
//    const XICPoints xicPointsSum = turboXIC.extractPointsXIC(100.0, 100.17, 2, 5);
//    QCOMPARE(xicPointsSum.scanNumbersVsIntensity.size(), 4);
//    QCOMPARE(xicPointsSum.scanNumbersVsIntensity.last(), 201.);

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

void TurboXICTests::mzMLMemoryTest() {

    ERR_INIT

    const QString file = "/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX43_A.raw.mzML";

//    QVector<int> test(2000);
//    QVector<QVector<int>> fdsa;
//    for (int i = 0; i < 200000; i++) {
//
//        if (i % 1000 == 0) {
//            qDebug() << i;
//        }
//
//        int testC[test.size()];
//        test.reserve(2000);
//
//        for (int j = 0; j < 2000; j++) {
//            testC[j] = j;
//        }
//
//        Eigen::VectorXi vec = Eigen::Map<Eigen::VectorXi>(testC, test.size());
////        std::vector<int> vec3(testC, testC + test.size());
////        fdsa.push_back(test);
//    }

    MsReaderMzML readerMzMl;
    e = readerMzMl.openFile(file);

    QMap<MzTargetKey, QMap<int, ScanPoints*>> pss;
    e = readerMzMl.collateMS2MzTargetFrames(&pss);

    QElapsedTimer et;
    et.start();

    QMap<MzTargetKey, TurboXIC*> turbos;
    for (auto it = pss.begin(); it != pss.end(); it++) {

        QMap<int, ScanPoints*> &ps = it.value();

        QElapsedTimer etLocal;
        etLocal.start();
        auto *turboXic = new TurboXIC();
        e = turboXic->init(ps);
        turbos.insert(it.key(), turboXic);
        qDebug() << it.key() << etLocal.elapsed() << "mSec";

        for (ScanPoints *p : ps) {
            ScanPoints().swap(*p);
        }
        QMap<int, ScanPoints*>().swap(ps);
    }

    readerMzMl.reset();
    qDebug() << et.elapsed();

    for (TurboXIC* ti : turbos) {
        delete ti;
    }

}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
