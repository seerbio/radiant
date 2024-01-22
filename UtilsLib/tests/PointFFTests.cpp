//
// Created by anichols on 11/07/2021.
//

#include "PointFF.h"

#include "MathUtils.h"

#include <QDebug>
#include <QtTest/QtTest>


class PointFFTests : public QObject
{
    Q_OBJECT

public:
    PointFFTests();
    ~PointFFTests() override = default;


private Q_SLOTS:

    static void xTest();
    static void yTest();
    static void rxTest();
    static void ryTest();
    static void equalOperatorTest();
    static void notEqualOperatorTest();

};


PointFFTests::PointFFTests() : QObject() {}

void PointFFTests::xTest() {

    PointFF pointFF;
    QVERIFY(MathUtils::tSame(pointFF.x(), 0.0f));

    PointFF pointDFLive(666.0f, 777.0f);
    QVERIFY(MathUtils::tSame(pointDFLive.x(), 666.0f));
}

void PointFFTests::yTest() {

    PointFF pointFF;
    QVERIFY(MathUtils::tSame(pointFF.y(), static_cast<float>(0.0)));

    PointFF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.y(), static_cast<float>(777.0)));

}

void PointFFTests::rxTest() {

    PointFF pointDFLive(666.0f, 777.0f);
    QVERIFY(MathUtils::tSame(pointDFLive.rx(), 666.0f));

    pointDFLive.rx() += 1.0;
    QVERIFY(MathUtils::tSame(pointDFLive.rx(), 667.0f));

}

void PointFFTests::ryTest() {

    PointFF pointDFLive(666.0, 777.0);
    QVERIFY(MathUtils::tSame(pointDFLive.ry(), static_cast<float>(777.0)));

    pointDFLive.ry() += 1.0;
    QVERIFY(MathUtils::tSame(pointDFLive.ry(), static_cast<float>(778.0)));
}

void PointFFTests::equalOperatorTest() {

    PointFF pointDF1(666.0, 777.0);
    PointFF pointDF2(667.0, 777.0);
    PointFF pointDF3(667.0, 778.0);

    QVERIFY(pointDF1 == pointDF1);
    QCOMPARE((pointDF1 == pointDF2), false);
    QCOMPARE((pointDF1 == pointDF3), false);

}

void PointFFTests::notEqualOperatorTest() {
    PointFF pointDF1(666.0, 777.0);
    PointFF pointDF2(667.0, 777.0);
    PointFF pointDF3(667.0, 778.0);

    QCOMPARE(pointDF1 != pointDF1, false);
    QVERIFY((pointDF1 != pointDF2));
    QVERIFY((pointDF1 != pointDF3));
}


QTEST_MAIN(PointFFTests)
#include "PointFFTests.moc"
