//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"

#include <QtTest/QtTest>

class ErrorUtilsTests : public QObject
{
    Q_OBJECT

public:
    ErrorUtilsTests();
    ~ErrorUtilsTests() override = default;

private Q_SLOTS:
    void toDoubleTest();
    void toIntTest();
    void emptyCheckReturnTest();
    void equalTest();
    void notEqualTest();
    void indexInRangeTest();
    void withinRangeTest();
    void aboveThresholdTest();
    void belowThresholdTest();
    void isTrueTest();
    void isFalseTest();
    void containsTest();
    void doesNotContainTest();

};


ErrorUtilsTests::ErrorUtilsTests()
{

}


void ErrorUtilsTests::toDoubleTest()
{
    ERR_INIT
    const QString valueToConvert = QStringLiteral("66.6");
    double value;
    e = ErrorUtils::toDouble(valueToConvert, &value, eValueError);
    QVERIFY(e == eNoError);

    const QString valueToConvertBad = QStringLiteral("66.6Da");
    e = ErrorUtils::toDouble(valueToConvertBad, &value, eValueError);
    QVERIFY(e == eValueError);

}


void ErrorUtilsTests::toIntTest()
{
    ERR_INIT
    const QString valueToConvert = QStringLiteral("666");
    int value;
    e = ErrorUtils::toInt(valueToConvert, &value, eValueError);
    QVERIFY(e == eNoError);

    const QString valueToConvertBad = QStringLiteral("666Da");
    e = ErrorUtils::toInt(valueToConvertBad, &value, eValueError);
    QVERIFY(e == eValueError);
}


void ErrorUtilsTests::emptyCheckReturnTest()
{
    const QVector<int> emptyInside;
    const QVector<int> numbersOfTheBeast = {666, 666};
    const QString nullString;
    ERR_INIT

    e = ErrorUtils::isNotEmpty(numbersOfTheBeast);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isNotEmpty(emptyInside, eValueError);
    QVERIFY(e == eValueError);

    eee_absorb;
    e = ErrorUtils::isNotEmpty(nullString);
    QVERIFY(e == eEmptyContainerError);

}


void ErrorUtilsTests::equalTest()
{
    ERR_INIT
    e = ErrorUtils::isEqual(66.6, 66.6);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isEqual(true, true);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isEqual(66.6, 66.661, eValueError);
    QVERIFY(e == eValueError);
}


void ErrorUtilsTests::isTrueTest()
{
    ERR_INIT
    e = ErrorUtils::isTrue(true);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isTrue(false, eValueError);
    QVERIFY(e == eValueError);
}


void ErrorUtilsTests::isFalseTest()
{
    ERR_INIT
    e = ErrorUtils::isFalse(false);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isFalse(true, eValueError);
    QVERIFY(e == eValueError);
}


void ErrorUtilsTests::indexInRangeTest()
{
    const QVector<int> numbersOfTheBeast = { 666, 666 };
    ERR_INIT

    e = ErrorUtils::isIndexInRange(numbersOfTheBeast, 0);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isIndexInRange(numbersOfTheBeast, 1);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isIndexInRange(numbersOfTheBeast,  -1, eValueError);
    QVERIFY(e == eValueError);

    eee_absorb;
    e = ErrorUtils::isIndexInRange(numbersOfTheBeast, 2);
    QVERIFY(e == eError);
}


void ErrorUtilsTests::notEqualTest()
{
    ERR_INIT

    e = ErrorUtils::isNotEqual(true, false);
    QVERIFY(e == eNoError);
}


void ErrorUtilsTests::withinRangeTest()
{
    ERR_INIT

    e = ErrorUtils::isWithinRange(666.6, 665.6, 667.6, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isWithinRange(666.6, 666.6, 667.6, ErrorUtilsParam::IncludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isWithinRange(667.6, 666.6, 667.6, ErrorUtilsParam::IncludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isWithinRange(666.6, 666.6, 667.6, ErrorUtilsParam::ExcludeThreshold, eValueError);
    QVERIFY(e == eValueError);

    eee_absorb;
    e = ErrorUtils::isWithinRange(667.6, 666.6, 667.6, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eError);
}


void ErrorUtilsTests::aboveThresholdTest()
{
    ERR_INIT

    e = ErrorUtils::isAboveThreshold(667, 666, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isAboveThreshold(666.6, 666.6, ErrorUtilsParam::IncludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isAboveThreshold(666, 667, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eError);

    eee_absorb;
    e = ErrorUtils::isAboveThreshold(666.6, 666.6, ErrorUtilsParam::ExcludeThreshold, eValueError);
    QVERIFY(e == eValueError);
}


void ErrorUtilsTests::belowThresholdTest()
{
    ERR_INIT

    e = ErrorUtils::isBelowThreshold(666, 667, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isBelowThreshold(666.6, 666.6, ErrorUtilsParam::IncludeThreshold);
    QVERIFY(e == eNoError);

    e = ErrorUtils::isBelowThreshold(667, 666, ErrorUtilsParam::ExcludeThreshold);
    QVERIFY(e == eError);

    eee_absorb;
    e = ErrorUtils::isBelowThreshold(666.6, 666.6, ErrorUtilsParam::ExcludeThreshold, eValueError);
    QVERIFY(e == eValueError);
}

void ErrorUtilsTests::containsTest() {

    ERR_INIT

    QVector<int> testVector{1, 2, 3, 4, 5};

    e = ErrorUtils::contains(3, testVector);
    QCOMPARE(e, eNoError);

    e = ErrorUtils::contains(666, testVector);
    QCOMPARE(e, eError);

    QMap<int, double> testMap{{1, 1.0}, {2, 1.0}, {3,1.0}};

    e = ErrorUtils::contains(3, testMap);
    QCOMPARE(e, eNoError);

    e = ErrorUtils::contains(666, testMap);
    QCOMPARE(e, eError);

    QHash<int, double> testHash{{1, 1.0}, {2, 1.0}, {3,1.0}};

    e = ErrorUtils::contains(3, testHash);
    QCOMPARE(e, eNoError);

    e = ErrorUtils::contains(666, testHash);
    QCOMPARE(e, eError);

}

void ErrorUtilsTests::doesNotContainTest() {

    ERR_INIT

    QVector<int> testVector{1, 2, 3, 4, 5};

    e = ErrorUtils::doesNotContain(3, testVector);
    QCOMPARE(e, eError);

    e = ErrorUtils::doesNotContain(666, testVector);
    QCOMPARE(e, eNoError);

    QMap<int, double> testMap{{1, 1.0}, {2, 1.0}, {3,1.0}};

    e = ErrorUtils::doesNotContain(3, testMap);
    QCOMPARE(e, eError);

    e = ErrorUtils::doesNotContain(666, testMap);
    QCOMPARE(e, eNoError);

    QHash<int, double> testHash{{1, 1.0}, {2, 1.0}, {3,1.0}};

    e = ErrorUtils::doesNotContain(3, testHash);
    QCOMPARE(e, eError);

    e = ErrorUtils::doesNotContain(666, testHash);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(ErrorUtilsTests)
#include "ErrorUtilsTests.moc"
