//
// Created by anichols on 11/07/2021.
//

#include "ObjectCSVWriters.h"

#include <QDebug>
#include <QtTest/QtTest>


class ObjectCSVWritersTests : public QObject
{
    Q_OBJECT

public:
    ObjectCSVWritersTests();
    ~ObjectCSVWritersTests() override = default;


private Q_SLOTS:

    static void writeReadTestVec();
    static void writeReadTestPairs();

};

ObjectCSVWritersTests::ObjectCSVWritersTests() : QObject() {}

void ObjectCSVWritersTests::writeReadTestVec() {

    ERR_INIT

    const QString writePath = QStringLiteral("writeMeAmadeus.csv");

    const QVector<int> vecToWrite = {6, 66, 666};

    e = ObjectCSVWriters::writeVectorToFile(vecToWrite, writePath);
    QCOMPARE(e, eNoError);

    QVector<int> vecToRead;
    e = ObjectCSVWriters::readVectorFromFile(writePath, &vecToRead);
    QCOMPARE(e, eNoError);
    QCOMPARE(vecToRead, vecToWrite);

}

void ObjectCSVWritersTests::writeReadTestPairs() {

    ERR_INIT

    const QString writePath = QStringLiteral("writeMeAmadeus2.csv");

    const QVector<QPair<int, float>> vecToWrite = {{6, 6.6}, {66, 66.6}, {666, 666.6}};

    e = ObjectCSVWriters::writeVectorToFile(vecToWrite, writePath);
    QCOMPARE(e, eNoError);

    QVector<QPair<int, float>> vecToRead;
    e = ObjectCSVWriters::readVectorFromFile(writePath, &vecToRead);
    QCOMPARE(e, eNoError);

    QCOMPARE(vecToRead, vecToWrite);
}


QTEST_MAIN(ObjectCSVWritersTests)
#include "ObjectCSVWritersTests.moc"
