//
// Created by anichols on 12/10/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILLBUILDER_H
#define PYTHIACPP_FEATUREFINDERHILLBUILDER_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHill.h"
#include "GlobalSettings.h"

#include <QScopedPointer>


using namespace Error;


struct ALGORITHMSLIB_EXPORTS FeatureFinderParameters {

    double tolerancePPM = -1.0;
    int skipScanCount = -1;
    int minScanCount = -1;

//    int filterLength = 7;
//    int smoothCount = 3;
//    double sigma = 1.0;
//    double signalToNoiseRatio = 3;


public:
    bool isValid() const {
        return tolerancePPM > 0.0
            && skipScanCount >= 0
            && minScanCount > 0;
//            && filterLength > 4
//            && smoothCount >0
//            && sigma > 0
//            && signalToNoiseRatio > 0;
    }

};


class ALGORITHMSLIB_EXPORTS FeatureFinderHillBuilder {


public:

    friend class FeatureFinderHillBuilderTests;

    FeatureFinderHillBuilder();
    ~FeatureFinderHillBuilder();

    Err init(const FeatureFinderParameters &featureFinderParameters);

    Err buildHills(
            const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber,
            QVector<FeatureFinderHill> *featureFinderHills
    );

    void setRunParallel(bool runParallel);

    static void writeHillsToBatmassMzMrtFile(
            const QMap<ScanNumber, double> &retentionTimeByScanNumber,
            const QVector<FeatureFinderHill> &featureFinderHills,
            const QString &destinationFilePath
    );


private:

    Err buildScanPointGroupsTest(
            const QVector<ScanPoints> &allScanPoints,
            QVector<QVector<QVector<double>>> *groupedMzVals,
            QVector<QVector<QVector<double>>> *groupedIntensityVals
    );

    Err connectCentroidsInGroupedMzValsTest(
            const QVector<QVector<QVector<double>>> &groupedMzVals,
            double tolerancePPM,
            QVector<QVector<int>> *connectedCentroidsVecs
    );


private:

    Q_DISABLE_COPY(FeatureFinderHillBuilder) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_FEATUREFINDERHILLBUILDER_H
