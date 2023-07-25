//
// Created by anichols on 7/1/23.
//

#ifndef PYTHIADIACPP_FRAGLIBIONRTREE_H
#define PYTHIADIACPP_FRAGLIBIONRTREE_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"

#include <QMap>

using namespace Error;

class MS2IonsSeparated;

struct FragLibIon {

    PeptideId peptideId = -1;
    double mz = -1.0;
    double rt = -1.0;
    IRT iRT = -1.0;
    double iMobility = -1.0;
    double intensity = -1.0;
    IonIndex ionIndex = -1;
    IonType ionType;
    Charge charge = -1;
    double frequencyPercent = -1.0;

    double mzSearched = -1.0;
    double intensitySearched = -1.0;
};


class ALGORITHMSLIB_EXPORTS FragLibIonRTree {

public:

    FragLibIonRTree();
    ~FragLibIonRTree();

    Err init(const QMap<PeptideStringWithMods, MS2IonsSeparated> &fragPreds);

    Err buildMzHashedVsFragLibIonFrequencePercentages(
            double ppmTolerance,
            QMap<MzHashed, FrequencyPercent> *mzHashVsFreqPct
            );

    Err addFrequencyPercentagesToFragLibIons(const QMap<MzHashed, FrequencyPercent> &mzHashVsFreqPct);

    Err updateFragLibIonsRTValues(const QString &iRTRecalibrationFilePath);

    [[nodiscard]] bool rtValsLoaded() const;

    Err getFragLibIons(
            double mzMin,
            double mzMax,
            double rtMin,
            double rtMax,
            QVector<FragLibIon> *foundFragLibIons
    );

    Err getFragLibIons(
            double mzMin,
            double mzMax,
            QVector<FragLibIon> *foundFragLibIons
    );

    Err getFragLibIonsByPeptideId(
            PeptideId peptideId,
            QVector<FragLibIon> *fragLibIons
            );

    Err getPeptideSequenceWithMods(
            PeptideId peptideId,
            PeptideStringWithMods *peptideStringWithMods
            );

private:

    Q_DISABLE_COPY(FragLibIonRTree) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_FRAGLIBIONRTREE_H
