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
    m_ionMobilityIndexVsDriftTime.clear();
    m_pointCount = 0;
    m_isInit = false;

    // Count only points that can actually enter the IM index.  The old code
    // reserved capacity for invalid-IM points even though it discarded them in
    // the population pass.
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

        m_pointCount += scanPoints.size();
        for (int pointIndex = 0; pointIndex < scanPoints.size(); ++pointIndex) {
            const float driftTime = timsbukIonMobilityOf(*alignedPointData, pointIndex);
            if (driftTime <= 0.0f) {
                continue;
            }

            ++mzBinVsCounts[mzBin(scanPoints.at(pointIndex).x())];
            m_ionMobilityIndexVsDriftTime.insert(
                ionMobilityIndexForDriftTime(driftTime),
                driftTime
                );
        }
    }

    m_mzBinVsPoints.reserve(mzBinVsCounts.size());
    for (auto countIt = mzBinVsCounts.constBegin(); countIt != mzBinVsCounts.constEnd(); ++countIt) {
        m_mzBinVsPoints[countIt.key()].reserve(countIt.value());
    }

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

        const FrameIndex frameIndex = msFrame.frameIndexFromScanNumber(alignedIt.key());
        for (int pointIndex = 0; pointIndex < scanPoints.size(); ++pointIndex) {
            const float driftTime = timsbukIonMobilityOf(*alignedPointData, pointIndex);
            if (driftTime <= 0.0f) {
                continue;
            }

            const ScanPoint &scanPoint = scanPoints.at(pointIndex);
            IndexedPoint indexedPoint;
            indexedPoint.mz = scanPoint.x();
            indexedPoint.intensity = scanPoint.y();
            indexedPoint.driftTime = driftTime;
            indexedPoint.frameIndex = frameIndex;
            indexedPoint.ionMobilityIndex = ionMobilityIndexForDriftTime(driftTime);
            m_mzBinVsPoints[mzBin(indexedPoint.mz)].push_back(indexedPoint);
        }
    }

    for (auto binIt = m_mzBinVsPoints.begin(); binIt != m_mzBinVsPoints.end(); ++binIt) {
        // MsFrame assigns a unique frame index to each scan in scan-number
        // order, and points are inserted in their original point-index order.
        // stable_sort therefore preserves the old final pointIndex tie-break
        // without storing pointIndex in every hot IndexedPoint.
        std::stable_sort(
            binIt.value().begin(),
            binIt.value().end(),
            [](const IndexedPoint &left, const IndexedPoint &right) {
                if (left.frameIndex != right.frameIndex) {
                    return left.frameIndex < right.frameIndex;
                }
                if (left.ionMobilityIndex != right.ionMobilityIndex) {
                    return left.ionMobilityIndex < right.ionMobilityIndex;
                }
                return left.mz < right.mz;
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

        const QVector<IndexedPoint> &points = binIt.value();
        auto beginIt = points.constBegin();
        auto endIt = points.constEnd();
        if (frameIndexMax > 0) {
            beginIt = std::lower_bound(
                points.constBegin(),
                points.constEnd(),
                frameIndexMin + 1,
                [](const IndexedPoint &point, FrameIndex frameIndex) {
                    return point.frameIndex < frameIndex;
                }
                );
            endIt = std::lower_bound(
                beginIt,
                points.constEnd(),
                frameIndexMax,
                [](const IndexedPoint &point, FrameIndex frameIndex) {
                    return point.frameIndex < frameIndex;
                }
                );
        }

        for (auto pointIt = beginIt; pointIt != endIt; ++pointIt) {
            const IndexedPoint &point = *pointIt;
            if (point.mz < mzMin || point.mz > mzMax) {
                continue;
            }
            if (!(ionMobilityMin <= point.driftTime && point.driftTime <= ionMobilityMax)) {
                continue;
            }

            XICPoint xicPoint;
            xicPoint.mz = point.mz;
            xicPoint.intensity = point.intensity;
            xicPoint.scanNumber = point.frameIndex;
            xicPoint.ionMobilityIndex = point.ionMobilityIndex;
            xicPoints.push_back(xicPoint);
        }
    }

    return xicPoints;
}
