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

    static void writeReadTest();

};

ObjectCSVWritersTests::ObjectCSVWritersTests() : QObject() {}

void ObjectCSVWritersTests::writeReadTest() {

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


QTEST_MAIN(ObjectCSVWritersTests)
#include "ObjectCSVWritersTests.moc"
