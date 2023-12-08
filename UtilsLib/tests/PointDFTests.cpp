//
// Created by anichols on 11/07/2021.
//

#include "PointDF.h"

#include "MathUtils.h"

#include <QDebug>
#include <QtTest/QtTest>


class PointDFTests : public QObject
{
    Q_OBJECT

public:
    PointDFTests();
    ~PointDFTests() override = default;


private Q_SLOTS:

    void xTest();
    void yTest();
    void rxTest();
    void ryTest();

};


PointDFTests::PointDFTests() : QObject() {}

void PointDFTests::xTest() {

    PointDF pointDF;
    QVERIFY(MathUtils::tSame(pointDF.x(), 0.0));

    PointDF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.x(), 666.0));
}

void PointDFTests::yTest() {

    PointDF pointDF;
    QVERIFY(MathUtils::tSame(pointDF.y(), static_cast<float>(0.0)));

    PointDF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.y(), static_cast<float>(777.0)));

}

void PointDFTests::rxTest() {

    PointDF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.rx(), 666.0));

    pointDFLive.rx() += 1.0;
    QVERIFY(MathUtils::tSame(pointDFLive.rx(), 667.0));

}

void PointDFTests::ryTest() {

    PointDF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.ry(), static_cast<float>(777.0)));

    pointDFLive.ry() += 1.0;
    QVERIFY(MathUtils::tSame(pointDFLive.ry(), static_cast<float>(778.0)));
}


QTEST_MAIN(PointDFTests)
#include "PointDFTests.moc"
