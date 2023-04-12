//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
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


private:

    QMap<int, QVector<QPointF>> buildPoints() const;


};


QMap<int, QVector<QPointF>> TurboXICTests::buildPoints() const {

    QMap<int, QVector<QPointF>> points;
    points.insert(1, {QPointF(100.11, 100.1), QPointF(200.11, 200.1)});
    points.insert(2, {QPointF(100.12, 100.2), QPointF(200.12, 200.2)});
    points.insert(3, {QPointF(100.13, 100.3), QPointF(200.13, 200.3)});
    points.insert(4, {QPointF(100.14, 100.4), QPointF(200.14, 200.4)});
    points.insert(5, {QPointF(100.15, 100.5), QPointF(100.151, 100.5), QPointF(200.15, 200.5)});


    return points;
}


void TurboXICTests::initTest() {

    const QMap<int, QVector<QPointF>> points = buildPoints();

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    e = turboXIC.init({});
    QCOMPARE(e, eError);

}


void TurboXICTests::extractPointsTest() {

    const QMap<int, QVector<QPointF>> points = buildPoints();

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    const XICPoints xicPoints = turboXIC.extractPoints(100.0, 100.13, 2, 4);
    QCOMPARE(xicPoints.size(), 2);
    QCOMPARE(xicPoints.firstKey(), 2);
    QCOMPARE(xicPoints.first(), 100.2);
    QCOMPARE(xicPoints.lastKey(), 3);
    QCOMPARE(xicPoints.last(), 100.3);

    const XICPoints xicPointsSum = turboXIC.extractPoints(100.0, 100.17, 2, 5);
    QCOMPARE(xicPointsSum.size(), 4);
    QCOMPARE(xicPointsSum.last(), 201.);

}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
