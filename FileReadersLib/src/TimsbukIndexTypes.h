//
// Created by Codex on 7/13/26.
//

#ifndef TIMSBUKINDEXTYPES_H
#define TIMSBUKINDEXTYPES_H

#include "FileReadersLib_Exports.h"

#include "GlobalSettings.h"

#include <QVector>
#include <QtGlobal>

using TimsbukCycleIndex = quint32;
using TimsbukWindowGroupId = int;

struct FILEREADERSLIB_EXPORTS TimsbukMs2WindowInfo {
    TimsbukWindowGroupId groupId = -1;
    float isolationMz = -1.0f;
    float isolationWidth = -1.0f;
    float collisionEnergy = -1.0f;
    QString relativePath;

    [[nodiscard]] bool isValid() const {
        return groupId >= 0 && isolationMz > 0.0f && isolationWidth > 0.0f;
    }
};

struct FILEREADERSLIB_EXPORTS TimsbukPeakRow {
    float mz = -1.0f;
    float intensity = -1.0f;
    float ionMobility = -1.0f;
    TimsbukCycleIndex cycleIndex = 0;

    [[nodiscard]] bool isValid() const {
        return mz > 0.0f && intensity >= 0.0f;
    }
};

struct FILEREADERSLIB_EXPORTS TimsbukLogicalScanDescriptor {
    ScanNumber scanNumber = -1;
    int msLevel = -1;
    TimsbukCycleIndex cycleIndex = 0;
    TimsbukWindowGroupId windowGroupId = -1;
    float scanTimeMilliseconds = -1.0f;
    float collisionEnergy = -1.0f;
    float precursorTargetMz = -1.0f;
    float isoWindowLower = -1.0f;
    float isoWindowUpper = -1.0f;
    MzTargetKey targetKey;
    int nativeFrameNumber = -1;

    [[nodiscard]] bool isMs1() const {
        return msLevel == 1;
    }

    [[nodiscard]] bool isMs2() const {
        return msLevel > 1;
    }
};

struct FILEREADERSLIB_EXPORTS TimsbukScanPointStore {
    ScanPoints scanPoints;
    // When populated, this stays index-aligned with scanPoints.
    QVector<float> ionMobilityByPoint;

    void clear() {
        scanPoints.clear();
        ionMobilityByPoint.clear();
    }

    [[nodiscard]] bool isEmpty() const {
        return scanPoints.isEmpty();
    }

    [[nodiscard]] bool hasIonMobility() const {
        return !ionMobilityByPoint.isEmpty();
    }

    [[nodiscard]] bool isAligned() const {
        return ionMobilityByPoint.isEmpty() || ionMobilityByPoint.size() == scanPoints.size();
    }

    [[nodiscard]] int pointCount() const {
        return scanPoints.size();
    }
};

struct FILEREADERSLIB_EXPORTS TimsbukLogicalScan {
    TimsbukLogicalScanDescriptor descriptor;
    TimsbukScanPointStore pointStore;

    void clear() {
        descriptor = {};
        pointStore.clear();
    }

    [[nodiscard]] bool isValid() const {
        return descriptor.scanNumber > 0 && pointStore.isAligned();
    }
};

[[nodiscard]] inline float timsbukMzOf(const ScanPoint &scanPoint) {
    return scanPoint.x();
}

[[nodiscard]] inline float timsbukIntensityOf(const ScanPoint &scanPoint) {
    return scanPoint.y();
}

[[nodiscard]] inline float timsbukIonMobilityOf(
    const TimsbukScanPointStore &pointStore,
    int pointIndex
    ) {

    if (pointIndex < 0 || pointIndex >= pointStore.ionMobilityByPoint.size()) {
        return -1.0f;
    }

    return pointStore.ionMobilityByPoint.at(pointIndex);
}

#endif // TIMSBUKINDEXTYPES_H
