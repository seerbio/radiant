//
// Created by anichols on 12/10/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILLBUILDER_H
#define PYTHIACPP_FEATUREFINDERHILLBUILDER_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHill.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

#include <QScopedPointer>


using namespace Error;


struct ALGORITHMSFFLIB_EXPORTS FeatureFinderParameters {

    //Hill building
    double tolerancePPM = -1.0;
    int skipScanCount = -1;
    int minScanCount = -1;

    //Hill Refinement Integration
    int filterLength = -1;
    int smoothCount = -1;
    double sigma = -1.0;
    double signalToNoiseRatio = -1.0;

    FeatureFinderParameters() = default;
    ~FeatureFinderParameters() = default;

    explicit FeatureFinderParameters(const PythiaParameters &pythiaParameters)
    : tolerancePPM(pythiaParameters.ms2ExtractionWidthPPM)
    , skipScanCount(pythiaParameters.skipScanCount)
    , minScanCount(pythiaParameters.minScanCount)
    , filterLength(pythiaParameters.filterLength)
    , smoothCount(pythiaParameters.smoothCount)
    , sigma(pythiaParameters.sigma)
    , signalToNoiseRatio(pythiaParameters.signalToNoiseRatio)
    {}

public:

    void printParams() const {

        qDebug() << "tolerancePPM" << tolerancePPM;
        qDebug() << "skipScanCount" << skipScanCount ;
        qDebug() << "minScanCount" << minScanCount;

        qDebug() << "filterLength" << filterLength;
        qDebug() << "smoothCount" << smoothCount;
        qDebug() << "Sigma" << sigma;
        qDebug() << "signalToNoiseRatio" << signalToNoiseRatio ;
    }

    [[nodiscard]] bool isValid() const {
        const bool isValid = tolerancePPM > 0.0
            && skipScanCount >= 0
            && minScanCount >= 1
            && filterLength >= 2
            && smoothCount > 0
            && sigma > 0.0
            && signalToNoiseRatio >= 1;

        if (!isValid) {
            printParams();
        }

        return isValid;
    }

};


class ALGORITHMSFFLIB_EXPORTS FeatureFinderHillBuilder {


public:

    friend class FeatureFinderHillBuilderTests;

    FeatureFinderHillBuilder();

    ~FeatureFinderHillBuilder();

    Err init(const FeatureFinderParameters &featureFinderParameters);

    Err buildHills(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    Err getHills(QVector<FeatureFinderHill*> *featureFinderHills);

    Err refineHills(bool useSmoothing);

    Err featureFinderHills(QVector<FeatureFinderHill> *featureFinderHills);

    void setRunParallel(bool runParallel);

    /*
     * destinationFilePath must end in .mzrt.csv extension
     */
    static Err writeHillsToBatmassMzMrtFile(
            const QMap<ScanNumber, double> &scanNumberVsScanTime,
            const QVector<FeatureFinderHill> &featureFinderHills,
            const QString &destinationFilePath
    );

private:

    Err buildScanPointGroupsTest(
            const QVector<ScanPoints*> &allScanPoints,
            QVector<QVector<QVector<double>>> *groupedMzVals,
            QVector<QVector<QVector<float>>> *groupedIntensityVals
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
