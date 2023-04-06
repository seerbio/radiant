//
// Created by anichols on 2/6/23.
//

#include "MsReaderParquet.h"

//#include "ErrorUtils.h"
//#include "StringUtils.h"
#include <QVector>
#include <QPointF>

#include <iostream>


bool MsReaderParquet::readFile(const std::string &fileURI) {


    return true;
}

bool MsReaderParquet::getScans(
        int msLevel,
        QMap<ScanNumber, ScanPoints> *scanPoints
        ) {

    if (m_parquetRows.empty()) {
        return false;
    }

    scanPoints->clear();

    for (const ParquetRow &pr : m_parquetRows) {

        if (pr.msLevel != msLevel) {
            continue;
        }

        (*scanPoints)[pr.scanNumber].push_back({pr.mz, pr.intensity});
    }

    return true;
}

bool MsReaderParquet::getScansInfo(QMap<ScanNumber, ScanInfo> *scansInfo) {

    if (m_parquetRows.empty()) {
        return false;
    }

    for (const ParquetRow &pr : m_parquetRows) {

        if (scansInfo->count(pr.scanNumber)) {
            continue;
        }

        ScanInfo scanInfo;
        scanInfo.scanNumber = pr.scanNumber;
        scanInfo.msLevel = pr.msLevel;
        scanInfo.retentionTime = pr.retentionTime;
        scanInfo.collisionEnergy = pr.collisionEnergy;
        scanInfo.mzTarget = pr.mzTarget;
        scanInfo.isoWindowLower = pr.isoWindowLower;
        scanInfo.isoWindowUpper = pr.isoWindowUpper;

        (*scansInfo)[scanInfo.scanNumber] = scanInfo;
    }

    return true;
}


bool MsReaderParquet::writeFile(const std::string &outputFilePath){


    return true;
}

void MsReaderParquet::clearParquetRows() {
    m_parquetRows.clear();
}

void MsReaderParquet::addParquetRowToParquetRows(const ParquetRow &parquetRow) {
    m_parquetRows.push_back(parquetRow);
}
