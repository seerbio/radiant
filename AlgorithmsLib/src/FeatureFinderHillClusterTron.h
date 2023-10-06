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

class MS2IonsSeparated;

struct HillsClustering {

    QVector<int> hillsMapIndexes;
    MonoOffset monoOffset = -1;
    double cosineSim = -1.0;
    Charge charge = -1;

};

struct HillsClusteringMS2 {
    double cosineSimSumWeighted = -1.0;
    FeatureFinderHillPlus apexFeatureFinderHillPlus;
    QVector<FeatureFinderHillPlus> correlatedHills;
};


class ALGORITHMSLIB_EXPORTS FeatureFinderHillClusterTron {

public:
    FeatureFinderHillClusterTron();
    ~FeatureFinderHillClusterTron() = default;

    Err init(
            const FeatureFinderParameters &params,
            int chargeMin,
            int chargeMax
            );

    Err init(const FeatureFinderParameters &params);

    Err clusterHillsMS1(
            const QVector<FeatureFinderHill> &featureFinderHills,
            bool isMs2Hills,
            QVector<HillsClustering> *hillClusterings,
            QMap<int, FeatureFinderHill> *featureFinderHillsMap
            );

    Err clusterHillsByBestMS2IonAnchor(
            const QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> &featureFinderHills,
            const QVector<MS2Ion> &ms2IonsAnchors,
            const MS2IonsSeparated &ms2IonsTandemPred,
            HillsClusteringMS2 *bestHillsClusteringMS2
    ) const;

    static Err writeClustersToMzRt(
            const QMap<ScanNumber, double> &scanNumberVsScanTime,
            const QVector<HillsClustering> &hillClustersByIndexs,
            const QMap<int, FeatureFinderHill> &featureFinderHillsMap,
            const QString &destinationFilePath
            );

    static Err calculateCosineSimBetweenHills(
            const FeatureFinderHill &h1,
            const FeatureFinderHill &h2,
            double *cosineSim
    );

private:

    FeatureFinderParameters m_params;

    int m_chargeMin;
    int m_chargeMax;

};


#endif //PYTHIADIACPP_FEATUREFINDERHILLCLUSTERTRON_H
