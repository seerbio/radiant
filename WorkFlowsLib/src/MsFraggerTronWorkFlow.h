//
// Created by anichols on 2/19/23.
//

#ifndef PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
#define PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibraryTron.h"
#include "FragmentLibraryRTree.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

class TandemScanIon;

using namespace Error;


struct PeptideIdIonFraggerResult {
    ScanNumber scanNumber = -1;
    PeptideId peptideId = -1;
    double searchedFragIonMz = -1.0;
    double foundTheoFragIonMz = -1.0;
    double intensity = -1.0;
    double ppmMzSearched = -1.0;
};

struct CHEMLIB_EXPORTS TallyPeptideId {
    ScanNumber scanNumber = -1;
    PeptideId peptideId = -1;
    Occurrence occurrence = 0;
    int intensityTotal = 0;
    double meanMzPPM = -1.0;
    QVector<double> scanIonMZs;
    QVector<double> theoFragMZs;
};

using FraggerKey = int;

class WORKFLOWSLIB_EXPORTS MsFraggerTronWorkFlow {

public:

    MsFraggerTronWorkFlow() = default;
    ~MsFraggerTronWorkFlow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &pepLibUri
            );

    Err processFile(const QString &mzmLFileURI);


private:

    Err buildRTrees(
            const QMap<int, QVector<TandemScanIon>> &tranchedTandemScanIons,
            QMap<int, FragmentLibraryRTree*> *rTreesByKey
            );

    Err fragScanIons(
            const QMap<FraggerKey, QVector<TandemScanIon>> &tranchedTandemScanIons,
            const QMap<FraggerKey, FragmentLibraryRTree*> &rTreesByKey,
            QHash<ScanNumber, QVector<PeptideIdIonFraggerResult>> *peptideIdIonFraggerResults
            );

    Err tallyPeptideIdsPerScan(
            const QHash<ScanNumber, QVector<PeptideIdIonFraggerResult>> &peptideIdIonFraggerResults,
            QMap<ScanNumber, QVector<TallyPeptideId>> *tallyItemsByScanNumber
            );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_pepLibUri;
};


#endif //PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
