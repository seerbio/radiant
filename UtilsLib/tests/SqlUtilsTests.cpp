//
// Created by anichols on 11/07/2021.
//

#include "SqlUtils.h"

#include <QDebug>
#include <QtTest/QtTest>


class SqlUtilsTests : public QObject
{
    Q_OBJECT

public:
    SqlUtilsTests();
    ~SqlUtilsTests() override = default;


private Q_SLOTS:

    void encodeBLOBTest();
    void decodeBLOBTest();

};


SqlUtilsTests::SqlUtilsTests() : QObject()
{
}


void SqlUtilsTests::encodeBLOBTest()
{

    const QVector<int> vec = {6, 6, 6, 666};
    const QByteArray vecEncode = SqlUtils::encodeBLOB(vec);
    const QString expectedResult = "\x06\x00\x00\x00\x06\x00\x00\x00\x06\x00\x00\x00\x9A\x02\x00\x00";
    QCOMPARE(vecEncode, expectedResult);
}


void SqlUtilsTests::decodeBLOBTest()
{

    const QVector<int> vec = {6, 6, 6, 666};
    const QByteArray vecEncode = SqlUtils::encodeBLOB(vec);
    const QString expectedResult = "\x06\x00\x00\x00\x06\x00\x00\x00\x06\x00\x00\x00\x9A\x02\x00\x00";
    QCOMPARE(vecEncode, expectedResult);

    const QVector<int> decodeVec = SqlUtils::decodeBLOB<int>(vecEncode);
    QCOMPARE(decodeVec, vec);

}


QTEST_MAIN(SqlUtilsTests)
#include "SqlUtilsTests.moc"
