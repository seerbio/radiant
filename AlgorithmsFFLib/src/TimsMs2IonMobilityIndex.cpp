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

    template<typename PointIterator>
    void frameRestrictedRange(
        PointIterator pointsBegin,
        PointIterator pointsEnd,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
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
            [](const auto &point, FrameIndex frameIndex) {
                return point.frameIndex < frameIndex;
            }
            );

        *endIt = std::lower_bound(
            *beginIt,
            pointsEnd,
            frameIndexMax,
            [](const auto &point, FrameIndex frameIndex) {
                return point.frameIndex < frameIndex;
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
    m_ionMobilityIndexVsDriftTime.clear();
    m_pointCount = 0;

    for (auto driftTimeIt = ionMobilityIndexVsDriftTime.constBegin();
         driftTimeIt != ionMobilityIndexVsDriftTime.constEnd();
         ++driftTimeIt) {
        m_ionMobilityIndexVsDriftTime.insert(
            driftTimeIt.key(),
            static_cast<float>(driftTimeIt.value())
            );
    }

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

            const float driftTime = driftTimeIt.value();
            const ScanPoints &scanPoints = mobilityIt.value();
            m_pointCount += scanPoints.size();
            for (const ScanPoint &scanPoint : scanPoints) {
                IndexedPoint indexedPoint;
                indexedPoint.mz = scanPoint.x();
                indexedPoint.intensity = scanPoint.y();
                indexedPoint.frameIndex = frameIndex;
                indexedPoint.ionMobilityIndex = mobilityIt.key();
                indexedPoint.driftTime = driftTime;
                m_mzBinVsPoints[mzBin(indexedPoint.mz)].push_back(indexedPoint);
            }
        }
    }

    for (auto binIt = m_mzBinVsPoints.begin(); binIt != m_mzBinVsPoints.end(); ++binIt) {
        std::sort(
            binIt.value().begin(),
            binIt.value().end(),
            [](const IndexedPoint &left, const IndexedPoint &right) {
                if (left.frameIndex != right.frameIndex) {
                    return left.frameIndex < right.frameIndex;
                }
                if (left.driftTime != right.driftTime) {
                    return left.driftTime < right.driftTime;
                }
                return left.mz < right.mz;
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

        const QVector<IndexedPoint> &points = binIt.value();
        QVector<IndexedPoint>::const_iterator beginIt;
        QVector<IndexedPoint>::const_iterator endIt;
        frameRestrictedRange(
            points.constBegin(),
            points.constEnd(),
            frameIndexMin,
            frameIndexMax,
            &beginIt,
            &endIt
            );

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

        const QVector<IndexedPoint> &points = binIt.value();
        QVector<IndexedPoint>::const_iterator beginIt;
        QVector<IndexedPoint>::const_iterator endIt;
        frameRestrictedRange(
            points.constBegin(),
            points.constEnd(),
            frameIndexMin,
            frameIndexMax,
            &beginIt,
            &endIt
            );

        for (auto pointIt = beginIt; pointIt != endIt; ++pointIt) {
            const IndexedPoint &point = *pointIt;
            if (point.mz < mzMin || point.mz > mzMax) {
                continue;
            }
            if (!(ionMobilityMin <= point.driftTime && point.driftTime <= ionMobilityMax)) {
                continue;
            }

            const double intensity = std::max(0.0f, point.intensity);
            (*mobilityProfile)[point.ionMobilityIndex] += intensity;

            if (point.intensity > *apexIntensity) {
                *apexIntensity = point.intensity;
                *apexDeltaAbs = std::abs(point.driftTime - ionMobilityCenter);
            }
        }
    }

    return !mobilityProfile->isEmpty();
}
