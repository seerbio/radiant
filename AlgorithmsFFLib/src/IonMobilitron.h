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
#include "XYMappermatic.h"

using namespace Error;

using IMPredicted = double;
using IMEmpirical = double;

class ALGORITHMSFFLIB_EXPORTS IonMobilitron {

public:

    IonMobilitron() = default;
    ~IonMobilitron() = default;

    Err init(const QVector<QPair<IMPredicted, IMEmpirical>> &imPredVsImEmpValuesSortedDiscScoreHiLo);

    Err predictIonMobilityIndex(
        float iIM,
        int *predictedIonMobilityIndex
        ) const;

    static Err assignIonMobilityValues(
        const PythiaParameters &pythiaParameters,
        const QVector<CandidateScores*> &candidateScorePairs,
        QMap<ScanNumber, FeatureFinderHillBuilder*> *scanNumberVsFeatureFinderHillBuildersPntrsTIMS
    );

private:

    XYMappermatic m_iIMtoIonMobilityIndexMapper;



};



#endif //IONMOBILITRON_H
