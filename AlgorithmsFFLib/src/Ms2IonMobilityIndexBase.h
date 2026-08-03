//
// Created by Codex on 7/29/26.
//

#ifndef PYTHIADIACPP_MS2IONMOBILITYINDEXBASE_H
#define PYTHIADIACPP_MS2IONMOBILITYINDEXBASE_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "TurboXIC.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS Ms2IonMobilityIndexBase {

public:

    virtual ~Ms2IonMobilityIndexBase() = default;

    [[nodiscard]] virtual bool isInit() const = 0;
    [[nodiscard]] virtual int pointCount() const = 0;
    [[nodiscard]] virtual bool driftTimeFromIonMobilityIndex(
        IonMobilityIndex ionMobilityIndex,
        float *driftTime
        ) const = 0;

    virtual XICPoints extractPointsXIC(
        float mzMin,
        float mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        float ionMobilityMin,
        float ionMobilityMax
        ) const = 0;
};

#endif //PYTHIADIACPP_MS2IONMOBILITYINDEXBASE_H
