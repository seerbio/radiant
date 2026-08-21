//
// Created by Codex on 7/23/26.
//

#ifndef PYTHIADIACPP_CENTROIDMS2IONMOBILITYINDEX_H
#define PYTHIADIACPP_CENTROIDMS2IONMOBILITYINDEX_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "Ms2IonMobilityIndexBase.h"
#include "TimsbukIndexTypes.h"

#include <QHash>
#include <QMap>
#include <QtGlobal>

class MsFrame;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS CentroidMs2IonMobilityIndex : public Ms2IonMobilityIndexBase {

public:

    Err init(
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        const QMap<ScanNumber, const TimsbukAlignedPointData*> &scanNumberVsAlignedPointData,
        const MsFrame &msFrame
        );

    [[nodiscard]] bool isInit() const override;
    [[nodiscard]] int pointCount() const override;
    [[nodiscard]] bool driftTimeFromIonMobilityIndex(
        IonMobilityIndex ionMobilityIndex,
        float *driftTime
        ) const override;

    XICPoints extractPointsXIC(
        float mzMin,
        float mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        float ionMobilityMin,
        float ionMobilityMax
        ) const override;

private:

    // Keep the values touched by extractPointsXIC() contiguous.  The previous
    // representation stored only indexes back into ScanPoints and aligned IM
    // vectors, which made every candidate query chase several unrelated
    // allocations for each point inspected.
    struct IndexedPoint {
        float mz = -1.0f;
        float intensity = -1.0f;
        float driftTime = -1.0f;
        FrameIndex frameIndex = -1;
        IonMobilityIndex ionMobilityIndex = -1;
    };

    QHash<int, QVector<IndexedPoint>> m_mzBinVsPoints;
    QMap<IonMobilityIndex, float> m_ionMobilityIndexVsDriftTime;
    int m_pointCount = 0;
    bool m_isInit = false;
};

#endif // PYTHIADIACPP_CENTROIDMS2IONMOBILITYINDEX_H
