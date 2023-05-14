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


private:

    PythiaParameters pythiaParameters() {

        PythiaParameters pythiaParameters;

        pythiaParameters.returnPSMTopN = 1;
        pythiaParameters.maxTandemPointCount = 2;
        pythiaParameters.ms2ExtractionWidthPPM = 12.0;
        pythiaParameters.precursorExtractionWindowThomsons = 1.0;
        pythiaParameters.chargeStateMin = 2;
        pythiaParameters.chargeStateMax = 3;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
        );

        return pythiaParameters;
    }
};

MsCalibratomaticTests::MsCalibratomaticTests() : QObject(){}

void MsCalibratomaticTests::execTests() {

    ERR_INIT

    const QString scoreVectorsFilePath
            = QDir(qApp->applicationDirPath()).filePath("test.prq");

    const int calPointK = 3;
    MsCalibratomatic calibratomatic;
    e = calibratomatic.init(
            pythiaParameters(),
            scoreVectorsFilePath,
            calPointK
            );
    QCOMPARE(e, eNoError);

    const QString mzMLFileURI
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(mzMLFileURI);
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanPoints> scanPoints = msReaderParquet.getScanPoints();

    QMap<ScanNumber, ScanPoints> scanPointsCal;
    e = calibratomatic.recalibratePoints(
            scanPoints,
            &scanPointsCal
            );
    QCOMPARE(e, eNoError);

    //TODO add more tests to actually test the corrected points.

}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
