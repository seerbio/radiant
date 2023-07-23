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

    LibraryBuilderWorkFlow();
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

    Err buildPeptideSequenceChargeKeyVsIsDecoy(const QVector<PeptidePredictionInput> &peptidePredictionInputs);

    Err buildPeptideStringWithModsVsIRT(const QVector<PeptidePredictionInput> &peptidePredictionInputs);

    Err writePredictionsToParquet(
            const QHash<PeptideSequenceChargeKey, TandemFragmentPredictotron::TandemPrediction> &tandemPredictionsAllCharges,
            const QString &outputFilePath
    );

private:

    QMap<Charge, QString> m_modelFilePaths;
    QMap<Charge, TandemFragmentPredictotron*> m_tandemPredictionModels;
    PythiaParameters m_pythiaParameters;
    QHash<PeptideSequenceChargeKey, bool> m_peptideSequenceChargeKeyVsIsDecoy;
    QHash<PeptideStringWithMods, IRT> m_peptideStringWithModsVsIRT;

    const int m_maxPeptideLength;

};


#endif //PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
