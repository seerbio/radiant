//
// Created by anichols on 5/26/23.
//

#ifndef PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H
#define PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHill.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"

using namespace Error;

struct HillsClustering {

    QVector<int> hillsMapIndexes;
    MonoOffset monoOffset = -1;
    double cosineSim = -1.0;
    Charge charge = -1;

};


class ALGORITHMSLIB_EXPORTS FeatureFinderHillClusterTron {

public:
    FeatureFinderHillClusterTron() = default;
    ~FeatureFinderHillClusterTron() = default;

    Err init(const FeatureFinderParameters &params);

    Err clusterHills(
            const QVector<FeatureFinderHill> &featureFinderHills,
            bool isMs2Hills,
            QVector<HillsClustering> *hillClusterings,
            QMap<int, FeatureFinderHill> *featureFinderHillsMap
            );


private:

    FeatureFinderParameters m_params;

};


#endif //PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H
