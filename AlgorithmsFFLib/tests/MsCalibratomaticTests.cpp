//
// Created by anichols on 11/07/2021.
//

#include "MsCalibratomatic.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>

class MsCalibratomaticTests : public QObject
{

Q_OBJECT

public:
    MsCalibratomaticTests();
    ~MsCalibratomaticTests() override = default;


private Q_SLOTS:

    void execTests();


};

MsCalibratomaticTests::MsCalibratomaticTests() : QObject(){}

void MsCalibratomaticTests::execTests() {
//    QSKIP("TODO: reenable with internal test data");

    ERR_INIT

    const QString reCalFile = QStringLiteral("/home/anichols/Desktop/Data/ConfigFiles/cal2.prq");

    MsCalibratomatic calibratomatic;
    e = calibratomatic.init(reCalFile);
    QCOMPARE(e, eNoError);

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = ParquetReader::read(reCalFile, &msCalibrationReaderRows);
    QCOMPARE(e, eNoError);

    e = calibratomatic.initMzOnly(msCalibrationReaderRows);
    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
