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
    double intensity = -1.0;
    double ppmMzSearched = -1.0;
};


class WORKFLOWSLIB_EXPORTS MsFraggerTronWorkFlow {

public:

    MsFraggerTronWorkFlow() = default;
    ~MsFraggerTronWorkFlow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri
            );

    Err processFile(const QString &mzmLFileURI);


private:

    Err buildRTrees(
            const QMap<int, QVector<TandemScanIon>> &tranchedTandemScanIons,
            QMap<int, FragmentLibraryRTree*> *rTreesByKey
            );

    Err fragScanIons(
            const QMap<int, QVector<TandemScanIon>> &tranchedTandemScanIons,
            const QMap<int, FragmentLibraryRTree*> &rTreesByKey,
            QMap<int, QVector<PeptideIdIonFraggerResult>> *peptideIdIonFraggerResults
            );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
};


#endif //PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
