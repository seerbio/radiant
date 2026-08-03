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

    struct SliceRef {
        const ScanPoints *scanPoints = nullptr;
        const TimsbukAlignedPointData *alignedPointData = nullptr;
        FrameIndex frameIndex = -1;
    };

    struct IndexedPointRef {
        quint32 sliceIndex = 0;
        quint32 pointIndex = 0;
        IonMobilityIndex ionMobilityIndex = -1;
    };

    [[nodiscard]] const SliceRef& sliceRefForPoint(const IndexedPointRef &pointRef) const;
    [[nodiscard]] const ScanPoint& scanPointForRef(const IndexedPointRef &pointRef) const;
    [[nodiscard]] float driftTimeForRef(const IndexedPointRef &pointRef) const;

    QHash<int, QVector<IndexedPointRef>> m_mzBinVsPoints;
    QVector<SliceRef> m_sliceRefs;
    QMap<IonMobilityIndex, float> m_ionMobilityIndexVsDriftTime;
    int m_pointCount = 0;
    bool m_isInit = false;
};

#endif // PYTHIADIACPP_CENTROIDMS2IONMOBILITYINDEX_H
