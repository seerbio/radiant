//
// Created by andrewnichols on 2/11/25.
//

#ifndef GENETICFEATURESELECTION_H
#define GENETICFEATURESELECTION_H


#include "AlgorithmsOptimizersLib_Exports.h"
#include "CandidateScores.h"
#include "GlobalSettings.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class ALGORITHMSOPTIMIZERSLIB_EXPORTS GeneticFeatureSelection {

public:
    GeneticFeatureSelection() = default;
    ~GeneticFeatureSelection() = default;

    static Err run(
        const PythiaParameters &pythiaParameters,
        const QVector<QString> &msDataFilePaths,
        const QVector<Features> &features,
        int maxGenerationIterations,
        int populationSize,
        double mutationRate,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        );

};



#endif //GENETICFEATURESELECTION_H
