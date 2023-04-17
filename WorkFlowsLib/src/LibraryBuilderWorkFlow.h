//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
#define PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "TandemFragmentPredictotron.h"


using namespace Error;


class Peptide;


class WORKFLOWSLIB_EXPORTS LibraryBuilderWorkFlow {

public:

    friend class LibraryBuilderWorkFlowTests;

    LibraryBuilderWorkFlow() = default;
    ~LibraryBuilderWorkFlow();

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &modelCharge1,
            const QString &modelCharge2,
            const QString &modelCharge3,
            const QString &modelCharge4
    );

    Err exec(
            const QString &peptidesCSVFilePath,
            QString *returnFilePath
            );


private:

    QMap<Charge, QString> m_modelFilePaths;
    QMap<Charge, TandemFragmentPredictotron*> m_tandemPredictionModels;
    PythiaParameters m_pythiaParameters;

};


#endif //PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
