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

    static void execTests();


};

MsCalibratomaticTests::MsCalibratomaticTests() : QObject(){}

void MsCalibratomaticTests::execTests() {

    ERR_INIT

    const QString reCalFile = QDir(qApp->applicationDirPath()).filePath("cal2.prq");

    MsCalibratomatic calibratomatic;

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = ParquetReader::read(reCalFile, &msCalibrationReaderRows);
    QCOMPARE(e, eNoError);

    e = calibratomatic.setMassCalibrationCoeffs(msCalibrationReaderRows, MSLevelEnum::MS2);
    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
