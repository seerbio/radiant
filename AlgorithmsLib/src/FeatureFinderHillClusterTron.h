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
    double cosineSimSum = -1.0;
    FeatureFinderHillPlus apexFeatureFinderHillPlus;
    QVector<FeatureFinderHillPlus> correlatedHills;
};

//TODO use this to deisotope feature finder hills.
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

    static Err clusterHillsByFrameIndex(
            const QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> &featureFinderHills,
            double mzMin,
            double mzMax,
            double cosineSimThreshold,
            HillsClusteringMS2 *bestHillsClusteringMS2
            );

private:

    FeatureFinderParameters m_params;

    double m_cosineSimThreshold;
    int m_chargeMin;
    int m_chargeMax;

};


#endif //PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H
