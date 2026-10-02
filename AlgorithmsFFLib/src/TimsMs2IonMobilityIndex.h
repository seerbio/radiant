//
// Created by Codex on 6/20/26.
//

#ifndef PYTHIADIACPP_TIMSMS2IONMOBILITYINDEX_H
#define PYTHIADIACPP_TIMSMS2IONMOBILITYINDEX_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "TurboXIC.h"

#include <QHash>
#include <QMap>

class MsFrame;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS TimsMs2IonMobilityIndex {

public:

    Err init(
        const QMap<FrameNumberTIMS, Ms2FrameTIMS> &frameNumberVsMs2FrameTims,
        const MsFrame &msFrame,
        const QMap<FrameIndex, double> &ionMobilityIndexVsDriftTime
        );

    [[nodiscard]] bool isInit() const;
    [[nodiscard]] int pointCount() const;
    [[nodiscard]] bool driftTimeFromIonMobilityIndex(
        IonMobilityIndex ionMobilityIndex,
        float *driftTime
        ) const;

    XICPoints extractPointsXIC(
        float mzMin,
        float mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        float ionMobilityMin,
        float ionMobilityMax
        ) const;

    bool extractMobilityProfile(
        float mzMin,
        float mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        float ionMobilityMin,
        float ionMobilityMax,
        float ionMobilityCenter,
        QMap<IonMobilityIndex, double> *mobilityProfile,
        float *apexIntensity,
        float *apexDeltaAbs
        ) const;

private:

    struct IndexedPoint {
        float mz = -1.0f;
        float intensity = 0.0f;
        FrameIndex frameIndex = -1;
        IonMobilityIndex ionMobilityIndex = -1;
        float driftTime = -1.0f;
    };

    QHash<int, QVector<IndexedPoint>> m_mzBinVsPoints;
    QMap<IonMobilityIndex, float> m_ionMobilityIndexVsDriftTime;
    int m_pointCount = 0;
    bool m_isInit = false;
};

#endif //PYTHIADIACPP_TIMSMS2IONMOBILITYINDEX_H
