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


struct HillsClusteringMS2 {

    QVector<int> hillsMapIndexes;
    double cosineSimSum = -1.0;
    Charge charge = -1;

};


class ALGORITHMSLIB_EXPORTS FeatureFinderHillClusterTron {

public:
    FeatureFinderHillClusterTron();
    ~FeatureFinderHillClusterTron() = default;

    Err init(
            const FeatureFinderParameters &params,
            double cosineSimThreshold,
            int chargeMin,
            int chargeMax
            );

    Err clusterHillsMS1(
            const QVector<FeatureFinderHill> &featureFinderHills,
            bool isMs2Hills,
            QVector<HillsClustering> *hillClusterings,
            QMap<int, FeatureFinderHill> *featureFinderHillsMap
            );

    static Err writeClustersToMzRt(
            const QMap<ScanNumber, double> &scanNumberVsScanTime,
            const QVector<HillsClustering> &hillClustersByIndexs,
            const QMap<int, FeatureFinderHill> &featureFinderHillsMap,
            const QString &destinationFilePath
            );

private:

    FeatureFinderParameters m_params;

    double m_cosineSimThreshold;
    int m_chargeMin;
    int m_chargeMax;

};


#endif //PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H
