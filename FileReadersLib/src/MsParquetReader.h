//
// Created by anichols on 2/6/23.
//

#ifndef SPARKDIA_PARQUETREADER_H
#define SPARKDIA_PARQUETREADER_H

//#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"

#include <QMap>
#include <vector>

//using namespace Error;


struct FILEREADERSLIB_EXPORTS ScanInfo {
    int msLevel = 1;
    int scanNumber = 1;
    double retentionTime = -1.0;
    double collisionEnergy = -1.0;
    double mzTarget = -1.0;
    double isoWindowLower = -1.0;
    double isoWindowUpper = -1.0;
};


struct FILEREADERSLIB_EXPORTS ParquetRow {

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


class FILEREADERSLIB_EXPORTS MsParquetReader {

public:

    MsParquetReader() = default;
    ~MsParquetReader() = default;

    static bool checkParquetStatus();

    bool readFile(const std::string &fileURI);

    bool writeFile(const std::string &outputFilePath);

    bool getScans(
            int msLevel,
            QMap<ScanNumber, ScanPoints> *scanPoints
            );

    bool getScansInfo(QMap<ScanNumber, ScanInfo> *scansInfo);

    void addParquetRowToParquetRows(const ParquetRow &parquetRow);

    void clearParquetRows();

private:

    std::vector<ParquetRow> m_parquetRows;

};


#endif //SPARKDIA_PARQUETREADER_H
