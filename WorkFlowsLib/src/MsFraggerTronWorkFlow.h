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

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
};


#endif //PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
