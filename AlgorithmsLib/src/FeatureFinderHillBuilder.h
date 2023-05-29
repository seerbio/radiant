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

    //Hill building
    double tolerancePPM = -1.0;
    int skipScanCount = -1;
    int minScanCount = -1;
    bool useMeanMz = false;

    //Hill Refinement Integration
    int filterLength = 5;
    int smoothCount = 1;
    double sigma = 1.0;
    double signalToNoiseRatio = 2;
    double cosineSimThreshold = 0.8;

    int chargeMin = 1;
    int chargeMax = 3;


public:

    bool isValid() const {
        const bool isValid = tolerancePPM > 0.0
            && skipScanCount >= 0
            && minScanCount >= 0
            && filterLength >= 2
            && smoothCount >0
            && sigma > 0
            && signalToNoiseRatio > 0;

        if (!isValid) {
           qDebug() << "tolerancePPM" << tolerancePPM;
           qDebug() << "skipScanCount" << skipScanCount ;
           qDebug() << "minScanCount" << minScanCount;
           qDebug() << "filterLength" << filterLength;
           qDebug() << "smoothCount" << smoothCount;
           qDebug() << "Sigma" << sigma;
           qDebug() << "signalToNoiseRatio" << signalToNoiseRatio ;
        }

        return isValid;
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

    /*
     * destinationFilePath must end in .mzrt.csv extension
     */
    static Err writeHillsToBatmassMzMrtFile(
            const QMap<ScanNumber, double> &scanNumberVsScanTime,
            const QVector<FeatureFinderHill> &featureFinderHills,
            const QString &destinationFilePath
    );

    Err refineHills(QVector<FeatureFinderHill> *featureFinderHills);

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
