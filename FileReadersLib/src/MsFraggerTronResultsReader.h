//
// Created by anichols on 2/22/23.
//

#ifndef PYTHIADIACPP_MSFRAGGERTRONRESULTSREADER_H
#define PYTHIADIACPP_MSFRAGGERTRONRESULTSREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"


using namespace Error;


struct RowToWrite {
    ScanNumber scanNumber = -1;
    PeptideId peptideId = -1;
    Occurrence occurrence = -1;
    int intensityTotal = 0;
    double meanMzPPM = -1.0;
    QString sequence;
    QString sequenceWithMods;
    QChar previousResidue;
    QChar postResidue;
    double mass = -1.0;
    bool isDecoy = false;
    QVector<double> scanIonMZs;
    QVector<double> theoFragIonMZs;
};


class FILEREADERSLIB_EXPORTS MsFraggerTronResultsReader {

public:

    static void writeToCsv(
            const QString &outputFilePath,
            QVector<RowToWrite> *rowsToWrite
            );

    static Err readCsv(
            const QString &firstPassCsvFilePath,
            QVector<RowToWrite> *rowsToWrite
            );


};


#endif //PYTHIADIACPP_MSFRAGGERTRONRESULTSREADER_H
