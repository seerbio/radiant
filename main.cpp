#include "src/ParquetReader.h"

#include "src/DeisotoperTandem.h"
#include "src/Error.h"
#include "src/ErrorUtils.h"
#include "src/StringUtils.h"

#include <iostream>
#include <map>

using namespace Error;

const std::string DOT_PARQUET = ".parquet";

struct CommandLineProcessotronVars {

    std::string fileUri;
    std::string outputUri;
    double ppmTol = -1.0;

};

Err replaceOriginalScansWithDeisotopedScans(
        const std::map<ScanNumber, ScanPoints> &deisotopedScanPoints,
        ParquetReader *parquetReader
        ) {

    ERR_INIT

    //NOTE: It is important to get scansInfo and MS1 Scans before clearing rows
    // in the next step.
    std::map<ScanNumber, ScanInfo> scanInfos;
    e = parquetReader->getScansInfo(&scanInfos); ree

    std::map<ScanNumber, ScanPoints> scanPointsMS1;
    e = parquetReader->getScans(1, &scanPointsMS1); ree;

    parquetReader->clearParquetRows();

    for (auto it = scanInfos.begin(); it != scanInfos.end(); it++) {

        const ScanNumber scanNumber = it->first;
        const ScanInfo &scanInfo = it->second;

        const ScanPoints &scanPoints = scanInfo.msLevel == 1
                ? scanPointsMS1.at(scanNumber)
                : deisotopedScanPoints.at(scanNumber);

        for (const ScanPoint &sp : scanPoints) {
            ParquetRow pr;
            pr.scanNumber = scanNumber;
            pr.mz = sp.x;
            pr.intensity = sp.y;
            pr.mzTarget = scanInfo.mzTarget;
            pr.msLevel = scanInfo.msLevel;
            pr.retentionTime = scanInfo.retentionTime;
            pr.collisionEnergy = scanInfo.collisionEnergy;
            pr.isoWindowLower = scanInfo.isoWindowLower;
            pr.isoWindowUpper = scanInfo.isoWindowUpper;

            parquetReader->addParquetRowToParquetRows(pr);
        }

    }

    ERR_RETURN
}

Err processCommandLine(
        int argc,
        char* argv[],
        CommandLineProcessotronVars *vars
) {

    ERR_INIT

    e = ErrorUtils::isTrue(argc > 1); ree;

    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];

        if (flag == "--file") {
            vars->fileUri = argv[i + 1];
        }

        else if (flag == "--outputUri") {
            vars->outputUri = argv[i + 1];
        }

        else if (flag == "--ppmTol") {
            vars->ppmTol = std::stod(argv[i + 1]);
        }
    }

    e = ErrorUtils::isTrue(vars->fileUri.find(DOT_PARQUET) != std::string::npos); ree;
    e = ErrorUtils::isNotEmpty(vars->fileUri); ree;
    e = ErrorUtils::isFalse(vars->ppmTol == -1); ree;

    ERR_RETURN
}

// ./SparkDIA --file /home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.parquet --ppmTol 25
int main(
        int argc,
        char* argv[]
        ) {

    ERR_INIT

    CommandLineProcessotronVars cliVars;

    e = processCommandLine(
            argc,
            argv,
            &cliVars);
    if (e != eNoError) {
        std::cout << "USAGE" << std::endl;
        std::cout << "./SparkDIA --file <path to parquet file> --ppmTol <ppm tolerance as an integer or decimal>" << std::endl;
        return 1;
    }

    std::cout << "Processing " << cliVars.fileUri << std::endl;

    const std::string &parquet_path = cliVars.fileUri;

    ParquetReader parquetReader;

    e = parquetReader.readFile(parquet_path);
    if (e != eNoError) {
        std::cout << errorMap.at(e) << std::endl;
        return 1;
    }

    std::map<ScanNumber, ScanPoints> scanPoints;
    e = parquetReader.getScans(2, &scanPoints);
    if (e != eNoError) {
        std::cout << errorMap.at(e) << std::endl;
        return 1;
    }

    std::map<ScanNumber, ScanPoints> deisotopedScanPoints;
    for (auto it = scanPoints.begin(); it != scanPoints.end(); it++) {

        const ScanNumber scanNumber = it->first;
        const ScanPoints scanPointsDeisotoped = DeisotoperTandem::deisotopeTandemScan(
                it->second,
                cliVars.ppmTol
                );

        deisotopedScanPoints.emplace(scanNumber, scanPointsDeisotoped);
    }

    e = replaceOriginalScansWithDeisotopedScans(
            deisotopedScanPoints,
            &parquetReader
            );
    if (e != eNoError) {
        std::cout << errorMap.at(e) << std::endl;
        return 1;
    }

    const std::string outputFilepath = cliVars.outputUri.empty()
            ? StringUtils::replaceSubstring(parquet_path,".parquet",".sansIsotopes.parquet")
            : cliVars.outputUri;

    e = parquetReader.writeFile(outputFilepath);
    if (e != eNoError) {
        std::cout << errorMap.at(e) << std::endl;
        return 1;
    }

    std::cout << "Processing complete, output located at: " << outputFilepath << std::endl;
    return 0;
}
