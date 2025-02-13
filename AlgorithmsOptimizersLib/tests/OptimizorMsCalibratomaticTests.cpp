#include <QtTest>
#include <QCoreApplication>

#include "OptimizorMsCalibratomatic.h"



#include <iostream>

class OptimizorMsCalibratomaticTests : public QObject {
Q_OBJECT

public:
    OptimizorMsCalibratomaticTests() = default;

    ~OptimizorMsCalibratomaticTests() override = default;

private slots:

    static void OptimizorMsCalibratomaticOptimize();

};


void OptimizorMsCalibratomaticTests::OptimizorMsCalibratomaticOptimize() {

    ERR_INIT

    const QString libFilePath = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib";
    const QString parametersFilePath = "/home/andrewnichols/Desktop/Data/ConfigFiles/test_params_decoys_v2_2_1_32cpu.pythiaConfig";
    const QVector<QString> msDataFilePaths = {
        "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML"
    };

    constexpr int maxGenerations = 50;
    e = OptimizorMsCalibratomatic::optimize(
        msDataFilePaths,
        libFilePath,
        parametersFilePath,
        maxGenerations
        );
    QCOMPARE(e, eNoError);

}

QTEST_MAIN(OptimizorMsCalibratomaticTests)

#include "OptimizorMsCalibratomaticTests.moc"
