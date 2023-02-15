//
// Created by anichols on 2/6/23.
//

#ifndef SPARKDIA_PARQUETREADER_H
#define SPARKDIA_PARQUETREADER_H

#include "ClassAliases.h"
#include "Error.h"

#include <vector>

using namespace Error;


struct ScanInfo {
    int msLevel = 1;
    int scanNumber = 1;
    double retentionTime = -1.0;
    double collisionEnergy = -1.0;
    double mzTarget = -1.0;
    double isoWindowLower = -1.0;
    double isoWindowUpper = -1.0;
};


struct ParquetRow {

    double mz = -1.0;
    float intensity = -1.0;
    int msLevel = -1;
    int scanNumber = 1;
    double retentionTime = -1.0;
    double collisionEnergy = -1.0;
    double mzTarget = -1.0;
    double isoWindowLower = -1.0;
    double isoWindowUpper = -1.0;
};


class ParquetReader {

public:

    ParquetReader() = default;
    ~ParquetReader() = default;

    static Err checkParquetStatus();

    Err readFile(const std::string &fileURI);

    Err writeFile(const std::string &outputFilePath);

    Err getScans(
            int msLevel,
            std::map<ScanNumber, ScanPoints> *scanPoints
            );

    Err getScansInfo(std::map<ScanNumber, ScanInfo> *scansInfo);

    void addParquetRowToParquetRows(const ParquetRow &parquetRow);

    void clearParquetRows();

private:

    std::vector<ParquetRow> m_parquetRows;

};


#endif //SPARKDIA_PARQUETREADER_H
