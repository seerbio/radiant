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

    const QString msDataFilePath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");

    const QString firstPassFilePath
            = QStringLiteral("/home/anichols/Repositories/PythiaDIACpp/FileReadersLib/tests/TestFiles/first_pass.csv");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");

    const int calPointK = 3;
    MsCalibratomatic calibratomatic;
    e = calibratomatic.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            firstPassFilePath,
            fragLibPath,
            msDataFilePath,
            calPointK
            );
    QCOMPARE(e, eNoError);

//    const QString mzMLFileURI
//            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");
//
//    MsReaderParquet msReaderParquet;
//    e = msReaderParquet.openFile(mzMLFileURI);
//    QCOMPARE(e, eNoError);
//
//    const QMap<ScanNumber, ScanPoints> scanPoints = msReaderParquet.getScanPoints();
//
//    QMap<ScanNumber, ScanPoints> scanPointsCal;
//    e = calibratomatic.recalibratePoints(
//            scanPoints,
//            &scanPointsCal
//            );
//    QCOMPARE(e, eNoError);

    //TODO add more tests to actually test the corrected points.

}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
