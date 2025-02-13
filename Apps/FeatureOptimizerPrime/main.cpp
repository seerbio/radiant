//
// Created by Drucifer on 12/31/2021.
//

#include "Error.h"

#include <QCoreApplication>
#include <QDirIterator>
#include <QElapsedTimer>

#include "OptimizorMsCalibratomatic.h"

using namespace Error;


int main(int argc, char *argv[]) {

    ERR_INIT

    const QString libFilePath = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib";
    const QString parametersFilePath = "/home/andrewnichols/Desktop/Data/ConfigFiles/test_params_decoys_v2_2_1_32cpu.pythiaConfig";
    const QVector<QString> msDataFilePaths = {
        // "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML",
        "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML"
    };

    constexpr int maxGenerations = 50;
    e = OptimizorMsCalibratomatic::optimize(
        msDataFilePaths,
        libFilePath,
        parametersFilePath,
        maxGenerations
        );
    if (e != eNoError) {
      qDebug() << "Error optimizing";
    }


    return 0;
}
