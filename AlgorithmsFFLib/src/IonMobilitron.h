//
// Created by andrewnichols on 1/7/25.
//

#ifndef IONMOBILITRON_H
#define IONMOBILITRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScores.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS IonMobilitron {

public:

    static Err assignIonMobilityValues(
        const PythiaParameters &pythiaParameters,
        const QVector<CandidateScores*> &candidateScorePairs,
        QMap<ScanNumber, FeatureFinderHillBuilder*> *scanNumberVsFeatureFinderHillBuildersPntrsTIMS
    );


private:



};



#endif //IONMOBILITRON_H
