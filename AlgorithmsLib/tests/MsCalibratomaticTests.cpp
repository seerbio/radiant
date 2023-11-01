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

    ERR_INIT

//    const QString reCalFile
//            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.prq.cached.pythiaCAL");
//
//    MsCalibratomatic calibratomatic;
//    e = calibratomatic.init(reCalFile);
//    QCOMPARE(e, eNoError);



    //TODO add more tests to actually test the corrected points.

}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
