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
#include <QtGlobal>

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

    // Points remain owned by the source Ms2FrameTIMS. The index stores
    // slice-local references so we do not duplicate mz/intensity payloads.
    struct SliceRef {
        const ScanPoints *scanPoints = nullptr;
        FrameIndex frameIndex = -1;
        IonMobilityIndex ionMobilityIndex = -1;
        float driftTime = -1.0f;
    };

    struct IndexedPointRef {
        quint32 sliceIndex = 0;
        quint32 pointIndex = 0;
    };

    [[nodiscard]] const SliceRef& sliceRefForPoint(const IndexedPointRef &pointRef) const;
    [[nodiscard]] const ScanPoint& scanPointForRef(const IndexedPointRef &pointRef) const;

    QHash<int, QVector<IndexedPointRef>> m_mzBinVsPoints;
    QVector<SliceRef> m_sliceRefs;
    QMap<IonMobilityIndex, float> m_ionMobilityIndexVsDriftTime;
    int m_pointCount = 0;
    bool m_isInit = false;
};

#endif //PYTHIADIACPP_TIMSMS2IONMOBILITYINDEX_H
