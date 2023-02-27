#include "ReCalibratomatic.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "MsFraggerTronResultsReader.h"

#include <QtTest>

#include <iostream>

#include <dlib/svm.h>

class ReCalibratomaticTests : public QObject
{
    Q_OBJECT

public:

    ReCalibratomaticTests() = default;
    ~ReCalibratomaticTests() override = default;

private slots:

    void svmTest();
    void recalibrateTandemScanIonsTests();

private:

    static QString fragFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("test" + S_GLOBAL_SETTINGS.DOT_FRAGLIB);
    }
};


void ReCalibratomaticTests::svmTest() {

    ReCalibratomatic::svmTest();

}

void ReCalibratomaticTests::recalibrateTandemScanIonsTests() {

    ERR_INIT

    const QString firstPassPSMs = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.psm.csv";

    QVector<RowToWrite> rowsToWrite;
    e = MsFraggerTronResultsReader::readCsv(
            firstPassPSMs,
            &rowsToWrite
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(rowsToWrite.size(), 1001);


    QVector<InputSVM> data;

    int counter = 0;
    for (const RowToWrite &rtw : rowsToWrite) {

        counter++;

        if (rtw.isDecoy) {
            break;
        }

        const QVector<double> &scanIonMZs = rtw.scanIonMZs;
        const QVector<double> &theoFragIonMZs = rtw.theoFragIonMZs;

        QCOMPARE(scanIonMZs.size(), theoFragIonMZs.size());

        for (int i = 0; i < scanIonMZs.size(); i++) {

            const double scanMz = scanIonMZs.at(i);
            const double theoMz = theoFragIonMZs.at(i);
            const double ppmDiff = 1e6 * (scanMz - theoMz) / theoMz;

            InputSVM is;
            is.scanNumber = rtw.scanNumber;
            is.mzScan = scanMz;
            is.mzTheo = theoMz;
            is.ppmDiff = ppmDiff;

            data.push_back(is);

        }

    }

    QCOMPARE(counter, 615);
    QCOMPARE(data.size(), 9714);

    ReCalibratomatic reCalibratomatic;
    e = reCalibratomatic.initSVM(data);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(ReCalibratomaticTests)

#include "ReCalibratomaticTests.moc"
