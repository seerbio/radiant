//
// Created by Codex on 7/23/26.
//

#include "CentroidMs2IonMobilityIndex.h"

#include "ErrorUtils.h"
#include "MsFrame.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr float CENTROID_MS2_MZ_BIN_WIDTH = 0.01f;
    constexpr float CENTROID_MS2_IM_INDEX_SCALE = 10000.0f;

    int mzBin(float mz) {
        return static_cast<int>(std::floor(mz / CENTROID_MS2_MZ_BIN_WIDTH));
    }

    IonMobilityIndex ionMobilityIndexForDriftTime(float driftTime) {
        return static_cast<IonMobilityIndex>(std::lround(driftTime * CENTROID_MS2_IM_INDEX_SCALE));
    }

    template<typename PointIterator, typename FrameIndexExtractor>
    void frameRestrictedRange(
        PointIterator pointsBegin,
        PointIterator pointsEnd,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        FrameIndexExtractor frameIndexForPoint,
        PointIterator *beginIt,
        PointIterator *endIt
        ) {

        *beginIt = pointsBegin;
        *endIt = pointsEnd;

        if (frameIndexMax <= 0) {
            return;
        }

        *beginIt = std::lower_bound(
            pointsBegin,
            pointsEnd,
            frameIndexMin + 1,
            [&](const auto &point, FrameIndex frameIndex) {
                return frameIndexForPoint(point) < frameIndex;
            }
            );

        *endIt = std::lower_bound(
            *beginIt,
            pointsEnd,
            frameIndexMax,
            [&](const auto &point, FrameIndex frameIndex) {
                return frameIndexForPoint(point) < frameIndex;
            }
            );
    }
}

Err CentroidMs2IonMobilityIndex::init(
    const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
    const QMap<ScanNumber, const TimsbukAlignedPointData*> &scanNumberVsAlignedPointData,
    const MsFrame &msFrame
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsAlignedPointData); ree;
    e = ErrorUtils::isTrue(msFrame.isValid()); ree;

    m_mzBinVsPoints.clear();
    m_sliceRefs.clear();
    m_ionMobilityIndexVsDriftTime.clear();
    m_pointCount = 0;
    m_isInit = false;

    m_sliceRefs.reserve(scanNumberVsAlignedPointData.size());

    QHash<int, int> mzBinVsCounts;
    for (auto alignedIt = scanNumberVsAlignedPointData.constBegin();
         alignedIt != scanNumberVsAlignedPointData.constEnd();
         ++alignedIt) {

        const auto scanPointIt = scanNumberVsScanPoints.constFind(alignedIt.key());
        if (scanPointIt == scanNumberVsScanPoints.constEnd()) {
            continue;
        }

        const TimsbukAlignedPointData *alignedPointData = alignedIt.value();
        if (alignedPointData == nullptr || !alignedPointData->isAlignedWith(scanPointIt.value())) {
            continue;
        }

        const ScanPoints &scanPoints = scanPointIt.value();
        if (scanPoints.isEmpty()) {
            continue;
        }

        SliceRef sliceRef;
        sliceRef.scanPoints = &scanPoints;
        sliceRef.alignedPointData = alignedPointData;
        sliceRef.frameIndex = msFrame.frameIndexFromScanNumber(alignedIt.key());
        m_sliceRefs.push_back(sliceRef);

        m_pointCount += scanPoints.size();
        for (int pointIndex = 0; pointIndex < scanPoints.size(); ++pointIndex) {
            ++mzBinVsCounts[mzBin(scanPoints.at(pointIndex).x())];

            const float driftTime = timsbukIonMobilityOf(*alignedPointData, pointIndex);
            if (driftTime <= 0.0f) {
                continue;
            }
            m_ionMobilityIndexVsDriftTime.insert(ionMobilityIndexForDriftTime(driftTime), driftTime);
        }
    }

    for (auto countIt = mzBinVsCounts.constBegin(); countIt != mzBinVsCounts.constEnd(); ++countIt) {
        m_mzBinVsPoints[countIt.key()].reserve(countIt.value());
    }

    for (quint32 sliceIndex = 0; sliceIndex < static_cast<quint32>(m_sliceRefs.size()); ++sliceIndex) {
        const SliceRef &sliceRef = m_sliceRefs.at(static_cast<int>(sliceIndex));
        for (quint32 pointIndex = 0; pointIndex < static_cast<quint32>(sliceRef.scanPoints->size()); ++pointIndex) {
            const float driftTime = timsbukIonMobilityOf(*sliceRef.alignedPointData, static_cast<int>(pointIndex));
            if (driftTime <= 0.0f) {
                continue;
            }

            IndexedPointRef pointRef;
            pointRef.sliceIndex = sliceIndex;
            pointRef.pointIndex = pointIndex;
            pointRef.ionMobilityIndex = ionMobilityIndexForDriftTime(driftTime);
            m_mzBinVsPoints[mzBin(sliceRef.scanPoints->at(static_cast<int>(pointIndex)).x())].push_back(pointRef);
        }
    }

    for (auto binIt = m_mzBinVsPoints.begin(); binIt != m_mzBinVsPoints.end(); ++binIt) {
        std::sort(
            binIt.value().begin(),
            binIt.value().end(),
            [this](const IndexedPointRef &left, const IndexedPointRef &right) {
                const SliceRef &leftSlice = sliceRefForPoint(left);
                const SliceRef &rightSlice = sliceRefForPoint(right);
                if (leftSlice.frameIndex != rightSlice.frameIndex) {
                    return leftSlice.frameIndex < rightSlice.frameIndex;
                }
                if (left.ionMobilityIndex != right.ionMobilityIndex) {
                    return left.ionMobilityIndex < right.ionMobilityIndex;
                }
                const float leftMz = scanPointForRef(left).x();
                const float rightMz = scanPointForRef(right).x();
                if (leftMz != rightMz) {
                    return leftMz < rightMz;
                }
                return left.pointIndex < right.pointIndex;
            }
            );
    }

    m_isInit = !m_mzBinVsPoints.isEmpty() && !m_ionMobilityIndexVsDriftTime.isEmpty();

    ERR_RETURN
}

bool CentroidMs2IonMobilityIndex::isInit() const {
    return m_isInit;
}

int CentroidMs2IonMobilityIndex::pointCount() const {
    return m_pointCount;
}

const CentroidMs2IonMobilityIndex::SliceRef& CentroidMs2IonMobilityIndex::sliceRefForPoint(
    const IndexedPointRef &pointRef
    ) const {

    return m_sliceRefs.at(static_cast<int>(pointRef.sliceIndex));
}

const ScanPoint& CentroidMs2IonMobilityIndex::scanPointForRef(
    const IndexedPointRef &pointRef
    ) const {

    const SliceRef &sliceRef = sliceRefForPoint(pointRef);
    return sliceRef.scanPoints->at(static_cast<int>(pointRef.pointIndex));
}

float CentroidMs2IonMobilityIndex::driftTimeForRef(const IndexedPointRef &pointRef) const {
    const SliceRef &sliceRef = sliceRefForPoint(pointRef);
    return timsbukIonMobilityOf(*sliceRef.alignedPointData, static_cast<int>(pointRef.pointIndex));
}

bool CentroidMs2IonMobilityIndex::driftTimeFromIonMobilityIndex(
    IonMobilityIndex ionMobilityIndex,
    float *driftTime
    ) const {

    if (driftTime == nullptr) {
        return false;
    }

    const auto it = m_ionMobilityIndexVsDriftTime.constFind(ionMobilityIndex);
    if (it == m_ionMobilityIndexVsDriftTime.constEnd()) {
        return false;
    }

    *driftTime = it.value();
    return true;
}

XICPoints CentroidMs2IonMobilityIndex::extractPointsXIC(
    float mzMin,
    float mzMax,
    FrameIndex frameIndexMin,
    FrameIndex frameIndexMax,
    float ionMobilityMin,
    float ionMobilityMax
    ) const {

    XICPoints xicPoints;

    if (!m_isInit || mzMax < mzMin || ionMobilityMax < ionMobilityMin) {
        return xicPoints;
    }

    const int binMin = mzBin(mzMin);
    const int binMax = mzBin(mzMax);
    for (int bin = binMin; bin <= binMax; ++bin) {
        const auto binIt = m_mzBinVsPoints.constFind(bin);
        if (binIt == m_mzBinVsPoints.constEnd()) {
            continue;
        }

        const QVector<IndexedPointRef> &points = binIt.value();
        QVector<IndexedPointRef>::const_iterator beginIt;
        QVector<IndexedPointRef>::const_iterator endIt;
        frameRestrictedRange(
            points.constBegin(),
            points.constEnd(),
            frameIndexMin,
            frameIndexMax,
            [this](const IndexedPointRef &pointRef) {
                return sliceRefForPoint(pointRef).frameIndex;
            },
            &beginIt,
            &endIt
            );

        for (auto pointIt = beginIt; pointIt != endIt; ++pointIt) {
            const IndexedPointRef &pointRef = *pointIt;
            const ScanPoint &scanPoint = scanPointForRef(pointRef);
            const float mz = scanPoint.x();
            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            const float driftTime = driftTimeForRef(pointRef);
            if (!(ionMobilityMin <= driftTime && driftTime <= ionMobilityMax)) {
                continue;
            }

            XICPoint xicPoint;
            xicPoint.mz = mz;
            xicPoint.intensity = scanPoint.y();
            xicPoint.scanNumber = sliceRefForPoint(pointRef).frameIndex;
            xicPoint.ionMobilityIndex = pointRef.ionMobilityIndex;
            xicPoints.push_back(xicPoint);
        }
    }

    return xicPoints;
}

bool CentroidMs2IonMobilityIndex::extractMobilityProfile(
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
    ) const {

    if (mobilityProfile == nullptr || apexIntensity == nullptr || apexDeltaAbs == nullptr) {
        return false;
    }

    mobilityProfile->clear();
    *apexIntensity = 0.0f;
    *apexDeltaAbs = std::max(std::abs(ionMobilityMax - ionMobilityCenter), std::abs(ionMobilityCenter - ionMobilityMin));

    if (!m_isInit || mzMax < mzMin || ionMobilityMax < ionMobilityMin) {
        return false;
    }

    const int binMin = mzBin(mzMin);
    const int binMax = mzBin(mzMax);
    for (int bin = binMin; bin <= binMax; ++bin) {
        const auto binIt = m_mzBinVsPoints.constFind(bin);
        if (binIt == m_mzBinVsPoints.constEnd()) {
            continue;
        }

        const QVector<IndexedPointRef> &points = binIt.value();
        QVector<IndexedPointRef>::const_iterator beginIt;
        QVector<IndexedPointRef>::const_iterator endIt;
        frameRestrictedRange(
            points.constBegin(),
            points.constEnd(),
            frameIndexMin,
            frameIndexMax,
            [this](const IndexedPointRef &pointRef) {
                return sliceRefForPoint(pointRef).frameIndex;
            },
            &beginIt,
            &endIt
            );

        for (auto pointIt = beginIt; pointIt != endIt; ++pointIt) {
            const IndexedPointRef &pointRef = *pointIt;
            const ScanPoint &scanPoint = scanPointForRef(pointRef);
            const float mz = scanPoint.x();
            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            const float driftTime = driftTimeForRef(pointRef);
            if (!(ionMobilityMin <= driftTime && driftTime <= ionMobilityMax)) {
                continue;
            }

            const float pointIntensity = scanPoint.y();
            const double intensity = std::max(0.0f, pointIntensity);
            (*mobilityProfile)[pointRef.ionMobilityIndex] += intensity;

            if (pointIntensity > *apexIntensity) {
                *apexIntensity = pointIntensity;
                *apexDeltaAbs = std::abs(driftTime - ionMobilityCenter);
            }
        }
    }

    return !mobilityProfile->isEmpty();
}
