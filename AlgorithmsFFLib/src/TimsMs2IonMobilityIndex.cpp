//
// Created by Codex on 6/20/26.
//

#include "TimsMs2IonMobilityIndex.h"

#include "ErrorUtils.h"
#include "MsFrame.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr float TIMS_MS2_MZ_BIN_WIDTH = 0.01f;

    int mzBin(float mz) {
        return static_cast<int>(std::floor(mz / TIMS_MS2_MZ_BIN_WIDTH));
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

Err TimsMs2IonMobilityIndex::init(
    const QMap<FrameNumberTIMS, Ms2FrameTIMS> &frameNumberVsMs2FrameTims,
    const MsFrame &msFrame,
    const QMap<FrameIndex, double> &ionMobilityIndexVsDriftTime
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(frameNumberVsMs2FrameTims); ree;
    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(ionMobilityIndexVsDriftTime); ree;

    m_mzBinVsPoints.clear();
    m_sliceRefs.clear();
    m_ionMobilityIndexVsDriftTime.clear();
    m_pointCount = 0;
    m_isInit = false;

    for (auto driftTimeIt = ionMobilityIndexVsDriftTime.constBegin();
         driftTimeIt != ionMobilityIndexVsDriftTime.constEnd();
         ++driftTimeIt) {
        m_ionMobilityIndexVsDriftTime.insert(
            driftTimeIt.key(),
            static_cast<float>(driftTimeIt.value())
            );
    }

    int sliceCount = 0;
    for (auto frameIt = frameNumberVsMs2FrameTims.constBegin();
         frameIt != frameNumberVsMs2FrameTims.constEnd();
         ++frameIt) {
        sliceCount += frameIt.value().size();
    }
    m_sliceRefs.reserve(sliceCount);

    QHash<int, int> mzBinVsCounts;
    for (auto frameIt = frameNumberVsMs2FrameTims.constBegin();
         frameIt != frameNumberVsMs2FrameTims.constEnd();
         ++frameIt) {

        const FrameIndex frameIndex = msFrame.frameIndexFromScanNumber(frameIt.key());

        for (auto mobilityIt = frameIt.value().constBegin();
             mobilityIt != frameIt.value().constEnd();
             ++mobilityIt) {

            const auto driftTimeIt = m_ionMobilityIndexVsDriftTime.constFind(mobilityIt.key());
            if (driftTimeIt == m_ionMobilityIndexVsDriftTime.constEnd()) {
                continue;
            }

            const ScanPoints &scanPoints = mobilityIt.value();
            if (scanPoints.isEmpty()) {
                continue;
            }

            SliceRef sliceRef;
            sliceRef.scanPoints = &scanPoints;
            sliceRef.frameIndex = frameIndex;
            sliceRef.ionMobilityIndex = mobilityIt.key();
            sliceRef.driftTime = driftTimeIt.value();
            m_sliceRefs.push_back(sliceRef);

            m_pointCount += scanPoints.size();
            for (const ScanPoint &scanPoint : scanPoints) {
                ++mzBinVsCounts[mzBin(scanPoint.x())];
            }
        }
    }

    for (auto countIt = mzBinVsCounts.constBegin(); countIt != mzBinVsCounts.constEnd(); ++countIt) {
        m_mzBinVsPoints[countIt.key()].reserve(countIt.value());
    }

    for (quint32 sliceIndex = 0; sliceIndex < static_cast<quint32>(m_sliceRefs.size()); ++sliceIndex) {
        const ScanPoints &scanPoints = *m_sliceRefs.at(static_cast<int>(sliceIndex)).scanPoints;
        for (quint32 pointIndex = 0; pointIndex < static_cast<quint32>(scanPoints.size()); ++pointIndex) {
            const ScanPoint &scanPoint = scanPoints.at(static_cast<int>(pointIndex));
            IndexedPointRef pointRef;
            pointRef.sliceIndex = sliceIndex;
            pointRef.pointIndex = pointIndex;
            m_mzBinVsPoints[mzBin(scanPoint.x())].push_back(pointRef);
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
                if (leftSlice.ionMobilityIndex != rightSlice.ionMobilityIndex) {
                    return leftSlice.ionMobilityIndex < rightSlice.ionMobilityIndex;
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

    m_isInit = m_pointCount > 0;

    ERR_RETURN
}

bool TimsMs2IonMobilityIndex::isInit() const {
    return m_isInit;
}

int TimsMs2IonMobilityIndex::pointCount() const {
    return m_pointCount;
}

const TimsMs2IonMobilityIndex::SliceRef& TimsMs2IonMobilityIndex::sliceRefForPoint(
    const IndexedPointRef &pointRef
    ) const {

    return m_sliceRefs.at(static_cast<int>(pointRef.sliceIndex));
}

const ScanPoint& TimsMs2IonMobilityIndex::scanPointForRef(
    const IndexedPointRef &pointRef
    ) const {

    const SliceRef &sliceRef = sliceRefForPoint(pointRef);
    return sliceRef.scanPoints->at(static_cast<int>(pointRef.pointIndex));
}

bool TimsMs2IonMobilityIndex::driftTimeFromIonMobilityIndex(
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

XICPoints TimsMs2IonMobilityIndex::extractPointsXIC(
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
            const SliceRef &sliceRef = sliceRefForPoint(pointRef);
            const ScanPoint &scanPoint = scanPointForRef(pointRef);
            const float mz = scanPoint.x();
            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            const float driftTime = sliceRef.driftTime;
            if (!(ionMobilityMin <= driftTime && driftTime <= ionMobilityMax)) {
                continue;
            }

            XICPoint xicPoint;
            xicPoint.mz = mz;
            xicPoint.intensity = scanPoint.y();
            xicPoint.scanNumber = sliceRef.frameIndex;
            xicPoint.ionMobilityIndex = sliceRef.ionMobilityIndex;
            xicPoints.push_back(xicPoint);
        }
    }

    return xicPoints;
}

bool TimsMs2IonMobilityIndex::extractMobilityProfile(
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
            const SliceRef &sliceRef = sliceRefForPoint(pointRef);
            const ScanPoint &scanPoint = scanPointForRef(pointRef);
            const float mz = scanPoint.x();
            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            const float driftTime = sliceRef.driftTime;
            if (!(ionMobilityMin <= driftTime && driftTime <= ionMobilityMax)) {
                continue;
            }

            const float pointIntensity = scanPoint.y();
            const double intensity = std::max(0.0f, pointIntensity);
            (*mobilityProfile)[sliceRef.ionMobilityIndex] += intensity;

            if (pointIntensity > *apexIntensity) {
                *apexIntensity = pointIntensity;
                *apexDeltaAbs = std::abs(driftTime - ionMobilityCenter);
            }
        }
    }

    return !mobilityProfile->isEmpty();
}
