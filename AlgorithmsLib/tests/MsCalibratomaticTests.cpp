//
// Created by anichols on 11/07/2021.
//

#include "MsCalibratomatic.h"
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
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores");

    const QString scansFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans");

    const QString fragLibFilePath
            = "/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib";

//    const QMap<QString, QString> filePaths = {
//            {scoreVectorsFilePath, scansFilePath}
//    };

    const QMap<QString, QString> filePaths = {
            {
                "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.595020.frameScores",
                "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.595020.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.615029.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.615029.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.404934.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.404934.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.414939.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.414939.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.424943.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.424943.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.434948.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.434948.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.444952.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.444952.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.454957.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.454957.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.464961.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.464961.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.484970.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.484970.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.494975.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.494975.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.504979.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.504979.frameScans"
            },
            {
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.514984.frameScores",
                    "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.514984.frameScans"
            },
    };


    const int calPointK = 25;
    MsCalibratomatic calibratomatic;
    e = calibratomatic.init(
            filePaths,
            pythiaParameters(),
            calPointK
            );
    QCOMPARE(e, eNoError);

    QVector<MsFrameScanPointRows> msFrameScanPointRows;
    ParquetReader::read(
            "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.595020.frameScans",
            &msFrameScanPointRows
    );
    QCOMPARE(e, eNoError);

    QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints;
    e = MsFrame::buildFrameIndexVsScanPoints(
            msFrameScanPointRows,
            &frameIndexVsScanPoints
            );
    QCOMPARE(e, eNoError);

    for (int i = 0; i < 62; i++) {
        qDebug() << i;
        QMap<FrameIndex, ScanPoints> recalFrameIndexVsScanPoints;
        e = calibratomatic.recalibratePoints(
                frameIndexVsScanPoints,
                &recalFrameIndexVsScanPoints
        );
        QCOMPARE(e, eNoError);

    }


    //TODO add more tests to actually test the corrected points.

}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
