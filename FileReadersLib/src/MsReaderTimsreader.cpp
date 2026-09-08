//
// Created by Codex on 9/8/26.
//

#include "MsReaderTimsreader.h"

#include "ErrorUtils.h"
#include "StringUtils.h"
#include "TimsreaderRunCatalog.h"

#include <QDebug>
#include <QDir>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#ifdef PYTHIA_HAVE_TIMSREADER
#include "timsreader.hpp"
#endif

namespace {

enum class TimsreaderTransformedColumn : int {
    ScanNumber = 0,
    MsLevel,
    ScanStartTimeSec,
    IsolationTargetMz,
    IsolationLowerOffset,
    IsolationUpperOffset,
    CollisionEnergy,
    IonMobilityWindowLower,
    IonMobilityWindowUpper,
    ScanLowerLimitMz,
    ScanUpperLimitMz,
    FrameIndex,
    CycleIndex,
    IsolationWindowId,
    IsolationRegionIndex,
    WindowGroupId,
    Masses,
    Intensities,
    IonMobilities,
    IonMobilityIndices,
};

constexpr int kFlattenedColumnCount = 18;
constexpr int kCentroidColumnCount = 20;

#ifdef PYTHIA_HAVE_TIMSREADER
constexpr size_t kLegacySidecarCentroidMaxPeaks = 20'000;
constexpr double kLegacySidecarCentroidMzTolerancePpm = 5.0;
constexpr double kLegacySidecarCentroidImPctTolerance = 3.0;
constexpr uint32_t kLegacySidecarCentroidEarlyStopIterations = 200;

tr_stream_config_t legacySidecarCentroidStreamConfig(const timsreader::Plan &plan) {
    tr_stream_config_t config = plan.default_stream_config();
    config.centroid_max_peaks = kLegacySidecarCentroidMaxPeaks;
    config.centroid_mz_tolerance_kind = TR_MZ_TOLERANCE_PPM;
    config.centroid_mz_tolerance_value = kLegacySidecarCentroidMzTolerancePpm;
    config.centroid_im_pct_tol = kLegacySidecarCentroidImPctTolerance;
    config.centroid_early_stop_iterations = kLegacySidecarCentroidEarlyStopIterations;
    config.centroid_window_cap_enabled = false;
    config.centroid_window_cap_max_peaks = 0;
    config.centroid_window_cap_window_da = 0.0f;
    return config;
}
#endif

int columnIndex(TimsreaderTransformedColumn column) {
    return static_cast<int>(column);
}

const ArrowArray *columnArray(const ArrowArray *recordBatch, TimsreaderTransformedColumn column) {
    if (recordBatch == nullptr
        || recordBatch->children == nullptr
        || columnIndex(column) >= recordBatch->n_children) {
        return nullptr;
    }

    return recordBatch->children[columnIndex(column)];
}

bool isValidAt(const ArrowArray *array, int64_t index) {
    if (array == nullptr) {
        return false;
    }
    if (array->null_count == 0 || array->buffers == nullptr || array->buffers[0] == nullptr) {
        return true;
    }

    const auto *validity = static_cast<const uint8_t*>(array->buffers[0]);
    const int64_t absoluteIndex = array->offset + index;
    return ((validity[absoluteIndex / 8] >> (absoluteIndex % 8)) & 0x1u) != 0;
}

template<typename T>
const T *valueBuffer(const ArrowArray *array) {
    if (array == nullptr || array->buffers == nullptr || array->n_buffers < 2) {
        return nullptr;
    }

    return static_cast<const T*>(array->buffers[1]);
}

template<typename T>
bool readNullablePrimitive(const ArrowArray *array, int64_t index, T *value) {
    if (value == nullptr || !isValidAt(array, index)) {
        return false;
    }

    const T *values = valueBuffer<T>(array);
    if (values == nullptr) {
        return false;
    }

    *value = values[array->offset + index];
    return true;
}

template<typename T>
T readRequiredPrimitive(const ArrowArray *array, int64_t index, T defaultValue = T()) {
    T value = defaultValue;
    readNullablePrimitive(array, index, &value);
    return value;
}

float halfBitsToFloat(uint16_t bits) {
    const uint32_t sign = static_cast<uint32_t>(bits & 0x8000u) << 16;
    uint32_t exponent = (bits >> 10) & 0x1fu;
    uint32_t mantissa = bits & 0x03ffu;
    uint32_t outBits = 0;

    if (exponent == 0) {
        if (mantissa == 0) {
            outBits = sign;
        }
        else {
            exponent = 127 - 15 + 1;
            while ((mantissa & 0x0400u) == 0) {
                mantissa <<= 1;
                --exponent;
            }
            mantissa &= 0x03ffu;
            outBits = sign | (exponent << 23) | (mantissa << 13);
        }
    }
    else if (exponent == 0x1fu) {
        outBits = sign | 0x7f800000u | (mantissa << 13);
    }
    else {
        outBits = sign | ((exponent + (127 - 15)) << 23) | (mantissa << 13);
    }

    float out = 0.0f;
    std::memcpy(&out, &outBits, sizeof(out));
    return out;
}

bool listRange(const ArrowArray *array, int64_t index, int32_t *begin, int32_t *end) {
    if (array == nullptr
        || begin == nullptr
        || end == nullptr
        || array->buffers == nullptr
        || array->n_buffers < 2
        || array->buffers[1] == nullptr
        || array->n_children < 1
        || array->children == nullptr
        || array->children[0] == nullptr
        || !isValidAt(array, index)) {
        return false;
    }

    const auto *offsets = static_cast<const int32_t*>(array->buffers[1]);
    const int64_t absoluteIndex = array->offset + index;
    *begin = offsets[absoluteIndex];
    *end = offsets[absoluteIndex + 1];
    return true;
}

void appendScanPointsFromLists(
    const ArrowArray *massesArray,
    const ArrowArray *intensitiesArray,
    int64_t rowIndex,
    ScanPoints *scanPoints
    ) {

    if (scanPoints == nullptr) {
        return;
    }

    scanPoints->clear();

    int32_t massBegin = 0;
    int32_t massEnd = 0;
    int32_t intensityBegin = 0;
    int32_t intensityEnd = 0;
    if (!listRange(massesArray, rowIndex, &massBegin, &massEnd)
        || !listRange(intensitiesArray, rowIndex, &intensityBegin, &intensityEnd)
        || (massEnd - massBegin) != (intensityEnd - intensityBegin)) {
        return;
    }

    const ArrowArray *massValuesArray = massesArray->children[0];
    const ArrowArray *intensityValuesArray = intensitiesArray->children[0];
    const float *massValues = valueBuffer<float>(massValuesArray);
    const float *intensityValues = valueBuffer<float>(intensityValuesArray);
    if (massValues == nullptr || intensityValues == nullptr) {
        return;
    }

    scanPoints->reserve(massEnd - massBegin);
    const int64_t massValueOffset = massValuesArray->offset;
    const int64_t intensityValueOffset = intensityValuesArray->offset;
    for (int32_t itemIndex = 0; itemIndex < massEnd - massBegin; ++itemIndex) {
        const float mz = massValues[massValueOffset + massBegin + itemIndex];
        const float intensity = intensityValues[intensityValueOffset + intensityBegin + itemIndex];
        scanPoints->push_back(ScanPoint(mz, intensity));
    }
}

void appendIonMobilitiesFromList(
    const ArrowArray *ionMobilitiesArray,
    int64_t rowIndex,
    QVector<float> *ionMobilities
    ) {

    if (ionMobilities == nullptr) {
        return;
    }

    ionMobilities->clear();

    int32_t begin = 0;
    int32_t end = 0;
    if (!listRange(ionMobilitiesArray, rowIndex, &begin, &end)) {
        return;
    }

    const ArrowArray *valuesArray = ionMobilitiesArray->children[0];
    const uint16_t *values = valueBuffer<uint16_t>(valuesArray);
    if (values == nullptr) {
        return;
    }

    ionMobilities->reserve(end - begin);
    const int64_t valueOffset = valuesArray->offset;
    for (int32_t itemIndex = begin; itemIndex < end; ++itemIndex) {
        ionMobilities->push_back(halfBitsToFloat(values[valueOffset + itemIndex]));
    }
}

float driftTimeFromWindowBounds(float lower, float upper) {
    if (lower <= 0.0f || upper <= 0.0f || upper < lower) {
        return -1.0f;
    }

    return (lower + upper) / 2.0f;
}

float minMz(const ScanPoints &scanPoints) {
    if (scanPoints.isEmpty()) {
        return std::numeric_limits<float>::max();
    }

    return std::min_element(
        scanPoints.begin(),
        scanPoints.end(),
        [](const ScanPoint &left, const ScanPoint &right) {
            return left.x() < right.x();
        }
        )->x();
}

float maxMz(const ScanPoints &scanPoints) {
    if (scanPoints.isEmpty()) {
        return -1.0f;
    }

    return std::max_element(
        scanPoints.begin(),
        scanPoints.end(),
        [](const ScanPoint &left, const ScanPoint &right) {
            return left.x() < right.x();
        }
        )->x();
}

}

class Q_DECL_HIDDEN MsReaderTimsreader::Private {
public:
    explicit Private(ImHandlingMode imHandlingMode)
        : m_imHandlingMode(imHandlingMode) {}

    void resetMaterializedState(MsReaderTimsreader *q) {
        runCatalog.close();
        alignedPointDataByScanNumber.clear();
        q->m_msScanInfo.clear();
        q->m_mzTargetVsScanInfosPntrs.clear();
        q->m_scanPoints.clear();
        q->m_scanNumberVsScanTime.clear();
        q->m_frameIndexVsDriftTime.clear();
        q->m_filePath.clear();
        q->m_mzMs1Min = std::numeric_limits<float>::max();
        q->m_mzMs1Max = -1.0f;
        q->m_mzMs2Min = std::numeric_limits<float>::max();
        q->m_mzMs2Max = -1.0f;
        q->setHasIonMobility(false);
#ifdef PYTHIA_HAVE_TIMSREADER
        ownedBatches.clear();
#endif
    }

    void rebuildIndexes(MsReaderTimsreader *q) {
        q->m_mzTargetVsScanInfosPntrs.clear();
        q->m_scanNumberVsScanTime.clear();
        q->m_mzMs1Min = std::numeric_limits<float>::max();
        q->m_mzMs1Max = -1.0f;
        q->m_mzMs2Min = std::numeric_limits<float>::max();
        q->m_mzMs2Max = -1.0f;

        for (auto it = q->m_msScanInfo.begin(); it != q->m_msScanInfo.end(); ++it) {
            MsScanInfo &msScanInfo = it.value();
            q->m_scanNumberVsScanTime.insert(msScanInfo.scanNumber, msScanInfo.scanTime);

            if (msScanInfo.msLevel > 1) {
                q->m_mzTargetVsScanInfosPntrs[msScanInfo.targetKey()].push_back(&msScanInfo);
            }

            const ScanPoints &scanPoints = q->m_scanPoints.value(msScanInfo.scanNumber);
            const float mzMin = minMz(scanPoints);
            const float mzMax = maxMz(scanPoints);
            if (msScanInfo.msLevel == 1) {
                q->m_mzMs1Min = std::min(q->m_mzMs1Min, mzMin);
                q->m_mzMs1Max = std::max(q->m_mzMs1Max, mzMax);
            }
            else {
                q->m_mzMs2Min = std::min(q->m_mzMs2Min, mzMin);
                q->m_mzMs2Max = std::max(q->m_mzMs2Max, mzMax);
            }
        }

        q->setHasIonMobility(m_imHandlingMode == ImHandlingMode::Centroid
            && !alignedPointDataByScanNumber.isEmpty());
    }

#ifdef PYTHIA_HAVE_TIMSREADER
    timsreader::Representation representation() const {
        switch (m_imHandlingMode) {
        case ImHandlingMode::Centroid:
            return timsreader::Representation::CentroidedFrameSpectrum;
        case ImHandlingMode::Summed:
            return timsreader::Representation::FlattenedWindowSpectrum;
        case ImHandlingMode::Raw4D:
            break;
        }

        return timsreader::Representation::DirectSubscans;
    }

    bool modeSupported() const {
        return m_imHandlingMode == ImHandlingMode::Centroid
            || m_imHandlingMode == ImHandlingMode::Summed;
    }

    Err ingestOwnedBatch(
        const timsreader::OwnedBatch &ownedBatch,
        MsReaderTimsreader *q
        ) {

        ERR_INIT

        const timsreader::BorrowedBatchView view = ownedBatch.export_view();
        const ArrowArray *recordBatch = view.array;
        e = ErrorUtils::isTrue(recordBatch != nullptr, eFileError); ree;

        const int expectedColumns = m_imHandlingMode == ImHandlingMode::Centroid
            ? kCentroidColumnCount
            : kFlattenedColumnCount;
        e = ErrorUtils::isTrue(recordBatch->n_children >= expectedColumns, eFileError); ree;

        const ArrowArray *scanNumberArray = columnArray(recordBatch, TimsreaderTransformedColumn::ScanNumber);
        const ArrowArray *msLevelArray = columnArray(recordBatch, TimsreaderTransformedColumn::MsLevel);
        const ArrowArray *scanTimeArray = columnArray(recordBatch, TimsreaderTransformedColumn::ScanStartTimeSec);
        const ArrowArray *isolationTargetMzArray = columnArray(recordBatch, TimsreaderTransformedColumn::IsolationTargetMz);
        const ArrowArray *isolationLowerOffsetArray = columnArray(recordBatch, TimsreaderTransformedColumn::IsolationLowerOffset);
        const ArrowArray *isolationUpperOffsetArray = columnArray(recordBatch, TimsreaderTransformedColumn::IsolationUpperOffset);
        const ArrowArray *collisionEnergyArray = columnArray(recordBatch, TimsreaderTransformedColumn::CollisionEnergy);
        const ArrowArray *ionMobilityWindowLowerArray = columnArray(recordBatch, TimsreaderTransformedColumn::IonMobilityWindowLower);
        const ArrowArray *ionMobilityWindowUpperArray = columnArray(recordBatch, TimsreaderTransformedColumn::IonMobilityWindowUpper);
        const ArrowArray *frameIndexArray = columnArray(recordBatch, TimsreaderTransformedColumn::FrameIndex);
        const ArrowArray *isolationRegionIndexArray = columnArray(recordBatch, TimsreaderTransformedColumn::IsolationRegionIndex);
        const ArrowArray *massesArray = columnArray(recordBatch, TimsreaderTransformedColumn::Masses);
        const ArrowArray *intensitiesArray = columnArray(recordBatch, TimsreaderTransformedColumn::Intensities);
        const ArrowArray *ionMobilitiesArray = m_imHandlingMode == ImHandlingMode::Centroid
            ? columnArray(recordBatch, TimsreaderTransformedColumn::IonMobilities)
            : nullptr;

        for (int64_t rowIndex = 0; rowIndex < recordBatch->length; ++rowIndex) {
            MsScanInfo msScanInfo;
            msScanInfo.scanNumber = static_cast<ScanNumber>(
                readRequiredPrimitive<uint32_t>(scanNumberArray, rowIndex) + 1u
                );
            msScanInfo.msLevel = static_cast<int>(
                readRequiredPrimitive<uint8_t>(msLevelArray, rowIndex)
                );
            msScanInfo.scanTime = static_cast<float>(
                readRequiredPrimitive<double>(scanTimeArray, rowIndex) / 60.0
                );

            readNullablePrimitive<float>(isolationTargetMzArray, rowIndex, &msScanInfo.precursorTargetMz);
            readNullablePrimitive<float>(isolationLowerOffsetArray, rowIndex, &msScanInfo.isoWindowLower);
            readNullablePrimitive<float>(isolationUpperOffsetArray, rowIndex, &msScanInfo.isoWindowUpper);
            readNullablePrimitive<float>(collisionEnergyArray, rowIndex, &msScanInfo.collisionEnergy);
            if (m_imHandlingMode == ImHandlingMode::Centroid) {
                readNullablePrimitive<float>(ionMobilityWindowLowerArray, rowIndex, &msScanInfo.ionMobilityWindowLower);
                readNullablePrimitive<float>(ionMobilityWindowUpperArray, rowIndex, &msScanInfo.ionMobilityWindowUpper);
                msScanInfo.ionMobilityDriftTime = driftTimeFromWindowBounds(
                    msScanInfo.ionMobilityWindowLower,
                    msScanInfo.ionMobilityWindowUpper
                    );
            }

            msScanInfo.nativeFrameNumber = static_cast<int>(
                readRequiredPrimitive<uint32_t>(frameIndexArray, rowIndex, 0u) + 1u
                );
            if (msScanInfo.msLevel > 1) {
                msScanInfo.nativeScanNumber = static_cast<int>(
                    readRequiredPrimitive<uint16_t>(isolationRegionIndexArray, rowIndex, 0u) + 1u
                    );
            }

            ScanPoints scanPoints;
            appendScanPointsFromLists(massesArray, intensitiesArray, rowIndex, &scanPoints);
            e = ErrorUtils::isNotEmpty(scanPoints, eFileError); ree;

            e = ErrorUtils::doesNotContain(msScanInfo.scanNumber, q->m_msScanInfo, eFileError); ree;
            e = ErrorUtils::doesNotContain(msScanInfo.scanNumber, q->m_scanPoints, eFileError); ree;
            q->m_msScanInfo.insert(msScanInfo.scanNumber, msScanInfo);
            q->m_scanPoints.insert(msScanInfo.scanNumber, scanPoints);

            if (m_imHandlingMode == ImHandlingMode::Centroid) {
                TimsbukAlignedPointData alignedPointData;
                appendIonMobilitiesFromList(
                    ionMobilitiesArray,
                    rowIndex,
                    &alignedPointData.ionMobilityByPoint
                    );
                if (alignedPointData.hasIonMobility()) {
                    e = ErrorUtils::isTrue(
                        alignedPointData.isAlignedWith(q->m_scanPoints.value(msScanInfo.scanNumber)),
                        eFileError
                        ); ree;
                    alignedPointDataByScanNumber.insert(msScanInfo.scanNumber, alignedPointData);
                }
            }
        }

        ERR_RETURN
    }

    Err loadPlan(
        const timsreader::Plan &plan,
        MsReaderTimsreader *q
        ) {

        ERR_INIT

        timsreader::Stream stream = m_imHandlingMode == ImHandlingMode::Centroid
            ? plan.open_stream(legacySidecarCentroidStreamConfig(plan))
            : plan.open_stream();
        while (true) {
            std::optional<timsreader::OwnedBatch> ownedBatch = stream.next_owned();
            if (!ownedBatch.has_value()) {
                break;
            }

            e = ingestOwnedBatch(*ownedBatch, q); ree;
            ownedBatches.push_back(std::move(*ownedBatch));
        }

        ERR_RETURN
    }

    Err materialize(const QString &filePath, MsReaderTimsreader *q) {
        ERR_INIT

        if (!modeSupported()) {
            rrr(eFunctionNotImplemented);
        }

        e = runCatalog.open(filePath); ree;

        try {
            timsreader::Run run(filePath.toStdString());
            e = loadPlan(
                run.make_plan(representation(), timsreader::Selector::ms1()),
                q
                ); ree;

            for (const TimsreaderIsolationWindowInfo &window : runCatalog.isolationWindows()) {
                e = loadPlan(
                    run.make_plan(
                        representation(),
                        timsreader::Selector::ms2_window(window.isolationWindowId)
                        ),
                    q
                    ); ree;
            }
        }
        catch (const std::exception &ex) {
            qDebug() << "timsreader failed to materialize Bruker run" << filePath << ex.what();
            resetMaterializedState(q);
            rrr(eFileError);
        }

        q->m_filePath = QDir::cleanPath(filePath);
        rebuildIndexes(q);

        ERR_RETURN
    }
#endif

    ImHandlingMode m_imHandlingMode = ImHandlingMode::Centroid;
    TimsreaderRunCatalog runCatalog;
    QMap<ScanNumber, TimsbukAlignedPointData> alignedPointDataByScanNumber;

#ifdef PYTHIA_HAVE_TIMSREADER
    std::vector<timsreader::OwnedBatch> ownedBatches;
#endif
};

MsReaderTimsreader::MsReaderTimsreader(ImHandlingMode imHandlingMode)
    : d_ptr(new Private(imHandlingMode)) {}

MsReaderTimsreader::~MsReaderTimsreader() = default;

Err MsReaderTimsreader::openFile(const QString &filePath) {
    ERR_INIT

    d_ptr->resetMaterializedState(this);
    e = ErrorUtils::fileExists(filePath); ree;

#ifdef PYTHIA_HAVE_TIMSREADER
    e = d_ptr->materialize(filePath, this); ree;
#else
    Q_UNUSED(filePath)
    rrr(eFunctionNotImplemented);
#endif

    ERR_RETURN
}

Err MsReaderTimsreader::openFile(
    const QString &filePath,
    const QString &columnToFilterBy,
    const QPair<double, double> &filterRange
    ) {

    ERR_INIT

    if (!columnToFilterBy.isEmpty()
        && !StringUtils::stringsMatch(columnToFilterBy, QStringLiteral("scanTime"), false)) {
        rrr(eFunctionNotImplemented);
    }

    e = openFile(filePath); ree;
    e = restrictScanTimeRange(
        static_cast<ScanTime>(filterRange.first),
        static_cast<ScanTime>(filterRange.second)
        ); ree;

    ERR_RETURN
}

Err MsReaderTimsreader::restrictScanTimeRange(ScanTime scanTimeMin, ScanTime scanTimeMax) {
    ERR_INIT

    e = MsReaderBase::restrictScanTimeRange(scanTimeMin, scanTimeMax); ree;

    for (auto it = d_ptr->alignedPointDataByScanNumber.begin();
         it != d_ptr->alignedPointDataByScanNumber.end();) {
        if (m_msScanInfo.contains(it.key())) {
            ++it;
        }
        else {
            it = d_ptr->alignedPointDataByScanNumber.erase(it);
        }
    }

    d_ptr->rebuildIndexes(this);

    ERR_RETURN
}

Err MsReaderTimsreader::getMzTargetScanPoints(
    const MzTargetKey &targetKey,
    QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
    ) {

    ERR_INIT

    scanNumberVsScanPoints->clear();

    e = ErrorUtils::isNotEmpty(targetKey); ree;
    e = ErrorUtils::isTrue(isInit()); ree;
    e = ErrorUtils::contains(targetKey, m_mzTargetVsScanInfosPntrs); ree;

    const QVector<MsScanInfo*> &targetMsScanInfos = m_mzTargetVsScanInfosPntrs.value(targetKey);
    e = ErrorUtils::isNotEmpty(targetMsScanInfos); ree;

    for (const MsScanInfo *msScanInfo : targetMsScanInfos) {
        e = ErrorUtils::isTrue(msScanInfo != nullptr, eFileError); ree;
        e = ErrorUtils::contains(msScanInfo->scanNumber, m_scanPoints); ree;
        scanNumberVsScanPoints->insert(msScanInfo->scanNumber, m_scanPoints.value(msScanInfo->scanNumber));
    }

    ERR_RETURN
}

Err MsReaderTimsreader::getMzTargetAlignedPointData(
    const MzTargetKey &targetKey,
    QMap<ScanNumber, const TimsbukAlignedPointData*> *scanNumberVsAlignedPointData
    ) const {

    ERR_INIT

    scanNumberVsAlignedPointData->clear();

    e = ErrorUtils::isNotEmpty(targetKey); ree;
    e = ErrorUtils::isTrue(isInit()); ree;
    e = ErrorUtils::contains(targetKey, m_mzTargetVsScanInfosPntrs); ree;

    const QVector<MsScanInfo*> &targetMsScanInfos = m_mzTargetVsScanInfosPntrs.value(targetKey);
    e = ErrorUtils::isNotEmpty(targetMsScanInfos); ree;

    for (const MsScanInfo *msScanInfo : targetMsScanInfos) {
        e = ErrorUtils::isTrue(msScanInfo != nullptr, eFileError); ree;
        const auto alignedPointDataIt
            = d_ptr->alignedPointDataByScanNumber.constFind(msScanInfo->scanNumber);
        if (alignedPointDataIt == d_ptr->alignedPointDataByScanNumber.constEnd()) {
            continue;
        }

        e = ErrorUtils::contains(msScanInfo->scanNumber, m_scanPoints); ree;
        e = ErrorUtils::isTrue(
            alignedPointDataIt.value().isAlignedWith(m_scanPoints.value(msScanInfo->scanNumber)),
            eFileError
            ); ree;
        scanNumberVsAlignedPointData->insert(msScanInfo->scanNumber, &alignedPointDataIt.value());
    }

    e = ErrorUtils::isNotEmpty(*scanNumberVsAlignedPointData); ree;

    ERR_RETURN
}

const TimsbukAlignedPointData *MsReaderTimsreader::alignedPointDataPntr(ScanNumber scanNumber) const {
    const auto alignedPointDataIt = d_ptr->alignedPointDataByScanNumber.constFind(scanNumber);
    if (alignedPointDataIt == d_ptr->alignedPointDataByScanNumber.constEnd()) {
        return nullptr;
    }

    return &alignedPointDataIt.value();
}

Err MsReaderTimsreader::closeFile() {
    d_ptr->resetMaterializedState(this);
    return eNoError;
}
