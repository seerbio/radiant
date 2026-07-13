//
// Created by Codex on 7/13/26.
//

#include <arrow/api.h>
#include <arrow/util/float16.h>
#include <parquet/arrow/reader.h>

// These need to come after Arrow/Parquet headers or Qt's `signals` macro breaks Arrow headers.
#include "MsReaderTimsbukIndex.h"

#include "ErrorUtils.h"
#include "StringUtils.h"

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QSet>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>


namespace {

    const QString TIMSBUK_INDEX_SUFFIX = QStringLiteral("idx");
    const QString TIMSBUK_METADATA_FILE_NAME = QStringLiteral("metadata.json");
    const QString TIMSBUK_BRUKER_INDEX_SUFFIX = QStringLiteral(".d.idx");
    const QString JSON_FIELD_BUCKET_SIZE = QStringLiteral("bucket_size");
    const QString JSON_FIELD_CREATED_AT = QStringLiteral("created_at");
    const QString JSON_FIELD_CYCLE_TO_RT_MS = QStringLiteral("cycle_to_rt_ms");
    const QString JSON_FIELD_GROUP_INFO = QStringLiteral("group_info");
    const QString JSON_FIELD_ID = QStringLiteral("id");
    const QString JSON_FIELD_MS1_PEAKS = QStringLiteral("ms1_peaks");
    const QString JSON_FIELD_MS2_WINDOW_GROUPS = QStringLiteral("ms2_window_groups");
    const QString JSON_FIELD_QUADRUPOLE_ISOLATION = QStringLiteral("quadrupole_isolation");
    const QString JSON_FIELD_RELATIVE_PATH = QStringLiteral("relative_path");
    const QString JSON_FIELD_VERSION = QStringLiteral("version");
    const QString JSON_FIELD_IM = QStringLiteral("im");
    const QString JSON_FIELD_MZ = QStringLiteral("mz");
    const float TIMSBUK_SCAN_TIME_MILLISECONDS_TO_MINUTES = 1.0f / 60000.0f;
    const float TIMSBUK_ION_MOBILITY_EPSILON = 0.0005f;

    struct TimsbukIsolationWindowSegment {
        TimsbukWindowGroupId groupId = -1;
        int windowIndex = -1;
        QString windowLabel;
        float precursorMzLower = -1.0f;
        float precursorMzUpper = -1.0f;
        float ionMobilityLower = -1.0f;
        float ionMobilityUpper = -1.0f;
        float collisionEnergy = -1.0f;
        QString relativePath;

        [[nodiscard]] bool isValid() const {
            return groupId >= 0
                && windowIndex >= 0
                && precursorMzLower > 0.0f
                && precursorMzUpper > precursorMzLower
                && ionMobilityUpper > ionMobilityLower;
        }

        [[nodiscard]] bool containsIonMobility(float ionMobility) const {
            return (ionMobilityLower - TIMSBUK_ION_MOBILITY_EPSILON) <= ionMobility
                && ionMobility <= (ionMobilityUpper + TIMSBUK_ION_MOBILITY_EPSILON);
        }

        [[nodiscard]] float precursorTargetMz() const {
            return (precursorMzLower + precursorMzUpper) / 2.0f;
        }

        [[nodiscard]] float isolationWindowLower() const {
            return precursorTargetMz() - precursorMzLower;
        }

        [[nodiscard]] float isolationWindowUpper() const {
            return precursorMzUpper - precursorTargetMz();
        }

        [[nodiscard]] float ionMobilityCenter() const {
            return (ionMobilityLower + ionMobilityUpper) / 2.0f;
        }

        [[nodiscard]] MzTargetKey targetKey() const {
            return MsScanInfo::targetKey(precursorMzLower, precursorMzUpper);
        }

        [[nodiscard]] int logicalWindowId() const {
            return (groupId * 100) + windowIndex;
        }
    };

    struct TimsbukPendingMs2Group {
        TimsbukMs2GroupMetadata metadata;
        QVector<TimsbukIsolationWindowSegment> isolationWindows;
        QVector<TimsbukScanPointStore> scanStores;

        [[nodiscard]] int cycleCount() const {
            return metadata.groupInfo.cycleToRtMilliseconds.size();
        }

        [[nodiscard]] int flatIndex(int windowIndex, TimsbukCycleIndex cycleIndex) const {
            return (windowIndex * cycleCount()) + static_cast<int>(cycleIndex);
        }

        [[nodiscard]] TimsbukScanPointStore *scanStore(int windowIndex, TimsbukCycleIndex cycleIndex) {
            return &scanStores[flatIndex(windowIndex, cycleIndex)];
        }

        [[nodiscard]] const TimsbukScanPointStore &scanStore(int windowIndex, TimsbukCycleIndex cycleIndex) const {
            return scanStores.at(flatIndex(windowIndex, cycleIndex));
        }
    };

    struct TimsbukPendingLogicalMs2Scan {
        TimsbukLogicalScan logicalScan;
        TimsbukIsolationWindowSegment isolationWindow;
    };

    QString cleanPath(const QString &filePath) {
        return QDir::cleanPath(filePath);
    }

    QString metadataFilePathForRoot(const QString &sidecarRootPath) {
        return QDir(sidecarRootPath).filePath(TIMSBUK_METADATA_FILE_NAME);
    }

    bool hasBrukerDirectorySuffix(const QFileInfo &fileInfo) {
        return fileInfo.isDir()
            && StringUtils::stringsMatch(
                fileInfo.suffix(),
                S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION,
                false
                );
    }

    bool hasIndexDirectorySuffix(const QFileInfo &fileInfo) {
        return fileInfo.isDir()
            && StringUtils::stringsMatch(
                fileInfo.suffix(),
                TIMSBUK_INDEX_SUFFIX,
                false
                );
    }

    bool directoryContainsMetadataFile(const QString &directoryPath) {
        const QFileInfo metadataFileInfo(metadataFilePathForRoot(directoryPath));
        return metadataFileInfo.exists() && metadataFileInfo.isFile();
    }

    QString inferSourceBrukerDirectoryPath(const QString &sidecarRootPath) {
        if (!sidecarRootPath.endsWith(TIMSBUK_BRUKER_INDEX_SUFFIX, Qt::CaseInsensitive)) {
            return {};
        }

        const QString brukerDirectoryPath = sidecarRootPath.left(sidecarRootPath.size() - 4);
        const QFileInfo brukerDirectoryInfo(brukerDirectoryPath);
        if (!hasBrukerDirectorySuffix(brukerDirectoryInfo)) {
            return {};
        }

        return cleanPath(brukerDirectoryPath);
    }

    Err validateSidecarRootPath(
        const QString &sidecarRootPath,
        const QString &contextPath
        ) {

        ERR_INIT

        const QFileInfo sidecarRootInfo(sidecarRootPath);
        if (!sidecarRootInfo.exists() || !sidecarRootInfo.isDir()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK sidecar root not found"
                     << "input_path" << contextPath
                     << "sidecar_root" << sidecarRootPath;
            rrr(eFileError);
        }

        const QString metadataFilePath = metadataFilePathForRoot(sidecarRootPath);
        const QFileInfo metadataFileInfo(metadataFilePath);
        if (!metadataFileInfo.exists() || !metadataFileInfo.isFile()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK sidecar root missing metadata.json"
                     << "input_path" << contextPath
                     << "sidecar_root" << sidecarRootPath
                     << "expected_metadata" << metadataFilePath;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parseRequiredStringField(
        const QJsonObject &jsonObject,
        const QString &fieldName,
        const QString &context,
        QString *value
        ) {

        ERR_INIT

        const QJsonValue jsonValue = jsonObject.value(fieldName);
        if (!jsonValue.isString()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK metadata field must be a string"
                     << "context" << context
                     << "field" << fieldName;
            rrr(eFileError);
        }

        *value = jsonValue.toString();
        e = ErrorUtils::isNotEmpty(*value, eFileError); ree;

        ERR_RETURN
    }

    Err parseRequiredNonNegativeIntField(
        const QJsonObject &jsonObject,
        const QString &fieldName,
        const QString &context,
        int *value
        ) {

        ERR_INIT

        const QJsonValue jsonValue = jsonObject.value(fieldName);
        if (!jsonValue.isDouble()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK metadata field must be numeric"
                     << "context" << context
                     << "field" << fieldName;
            rrr(eFileError);
        }

        const double number = jsonValue.toDouble();
        if (number < 0.0 || std::floor(number) != number) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK metadata field must be a non-negative integer"
                     << "context" << context
                     << "field" << fieldName
                     << "value" << number;
            rrr(eFileError);
        }

        *value = static_cast<int>(number);

        ERR_RETURN
    }

    Err parseCycleToRtMilliseconds(
        const QJsonValue &jsonValue,
        const QString &context,
        QVector<float> *cycleToRtMilliseconds
        ) {

        ERR_INIT

        cycleToRtMilliseconds->clear();

        if (jsonValue.isArray()) {
            const QJsonArray jsonArray = jsonValue.toArray();
            if (jsonArray.isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK cycle_to_rt_ms array is empty"
                         << "context" << context;
                rrr(eFileError);
            }

            cycleToRtMilliseconds->reserve(jsonArray.size());
            for (const QJsonValue &entryValue : jsonArray) {
                if (!entryValue.isDouble()) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK cycle_to_rt_ms entry must be numeric"
                             << "context" << context;
                    rrr(eFileError);
                }

                cycleToRtMilliseconds->push_back(static_cast<float>(entryValue.toDouble()));
            }

            ERR_RETURN
        }

        if (jsonValue.isObject()) {
            const QJsonObject jsonObject = jsonValue.toObject();
            if (jsonObject.isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK cycle_to_rt_ms object is empty"
                         << "context" << context;
                rrr(eFileError);
            }

            QMap<int, float> cycleIndexVsRtMilliseconds;
            int cycleIndexMax = -1;
            for (auto it = jsonObject.constBegin(); it != jsonObject.constEnd(); ++it) {
                bool ok = false;
                const int cycleIndex = it.key().toInt(&ok);
                if (!ok || cycleIndex < 0 || !it.value().isDouble()) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK cycle_to_rt_ms object entry invalid"
                             << "context" << context
                             << "key" << it.key();
                    rrr(eFileError);
                }

                cycleIndexVsRtMilliseconds.insert(
                    cycleIndex,
                    static_cast<float>(it.value().toDouble())
                    );
                cycleIndexMax = std::max(cycleIndexMax, cycleIndex);
            }

            cycleToRtMilliseconds->fill(-1.0f, cycleIndexMax + 1);
            for (auto it = cycleIndexVsRtMilliseconds.constBegin(); it != cycleIndexVsRtMilliseconds.constEnd(); ++it) {
                (*cycleToRtMilliseconds)[it.key()] = it.value();
            }

            const bool cycleIndicesAreContiguous = std::all_of(
                cycleToRtMilliseconds->begin(),
                cycleToRtMilliseconds->end(),
                [](float rtMilliseconds){ return rtMilliseconds >= 0.0f; }
                );
            if (!cycleIndicesAreContiguous) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK cycle_to_rt_ms object keys must be contiguous"
                         << "context" << context;
                rrr(eFileError);
            }

            ERR_RETURN
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMSBUK cycle_to_rt_ms must be an array or object"
                 << "context" << context;
        rrr(eFileError);
    }

    Err validateRelativeDataPath(
        const QString &sidecarRootPath,
        const QString &relativePath,
        const QString &context
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(relativePath, eFileError); ree;
        if (QDir::isAbsolutePath(relativePath)) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK metadata relative_path must be relative"
                     << "context" << context
                     << "relative_path" << relativePath;
            rrr(eFileError);
        }

        const QString resolvedPath = cleanPath(QDir(sidecarRootPath).filePath(relativePath));
        const QFileInfo fileInfo(resolvedPath);
        if (!fileInfo.exists() || !fileInfo.isFile()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK data file missing"
                     << "context" << context
                     << "relative_path" << relativePath
                     << "resolved_path" << resolvedPath;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parsePeakGroupMetadata(
        const QJsonObject &jsonObject,
        const QString &sidecarRootPath,
        const QString &context,
        TimsbukPeakGroupMetadata *peakGroupMetadata
        ) {

        ERR_INIT

        peakGroupMetadata->clear();

        e = parseRequiredStringField(
            jsonObject,
            JSON_FIELD_RELATIVE_PATH,
            context,
            &peakGroupMetadata->relativePath
            ); ree;
        e = validateRelativeDataPath(
            sidecarRootPath,
            peakGroupMetadata->relativePath,
            context
            ); ree;
        e = parseCycleToRtMilliseconds(
            jsonObject.value(JSON_FIELD_CYCLE_TO_RT_MS),
            context,
            &peakGroupMetadata->cycleToRtMilliseconds
            ); ree;
        e = parseRequiredNonNegativeIntField(
            jsonObject,
            JSON_FIELD_BUCKET_SIZE,
            context,
            &peakGroupMetadata->bucketSize
            ); ree;

        if (!peakGroupMetadata->isValid()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK peak group metadata invalid"
                     << "context" << context;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parseQuadrupoleIsolation(
        const QJsonValue &jsonValue,
        const QString &context,
        QJsonArray *quadrupoleIsolation
        ) {

        ERR_INIT

        *quadrupoleIsolation = QJsonArray();

        if (jsonValue.isArray()) {
            *quadrupoleIsolation = jsonValue.toArray();
        }
        else if (jsonValue.isObject()) {
            quadrupoleIsolation->push_back(jsonValue.toObject());
        }
        else {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK quadrupole_isolation must be an array or object"
                     << "context" << context;
            rrr(eFileError);
        }

        if (quadrupoleIsolation->isEmpty()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK quadrupole_isolation is empty"
                     << "context" << context;
            rrr(eFileError);
        }

        for (const QJsonValue &entryValue : *quadrupoleIsolation) {
            if (!entryValue.isObject() || entryValue.toObject().isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK quadrupole_isolation entry invalid"
                         << "context" << context;
                rrr(eFileError);
            }
        }

        ERR_RETURN
    }

    Err parseMs2GroupMetadata(
        const QJsonObject &jsonObject,
        const QString &sidecarRootPath,
        TimsbukMs2GroupMetadata *ms2GroupMetadata
        ) {

        ERR_INIT

        ms2GroupMetadata->clear();

        const QString context = QStringLiteral("ms2_window_groups[]");

        e = parseRequiredNonNegativeIntField(
            jsonObject,
            JSON_FIELD_ID,
            context,
            &ms2GroupMetadata->groupId
            ); ree;
        e = parseQuadrupoleIsolation(
            jsonObject.value(JSON_FIELD_QUADRUPOLE_ISOLATION),
            context,
            &ms2GroupMetadata->quadrupoleIsolation
            ); ree;

        const QJsonValue groupInfoValue = jsonObject.value(JSON_FIELD_GROUP_INFO);
        if (!groupInfoValue.isObject()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK group_info must be an object"
                     << "context" << context
                     << "group_id" << ms2GroupMetadata->groupId;
            rrr(eFileError);
        }

        e = parsePeakGroupMetadata(
            groupInfoValue.toObject(),
            sidecarRootPath,
            QStringLiteral("ms2_window_groups[%1].group_info").arg(ms2GroupMetadata->groupId),
            &ms2GroupMetadata->groupInfo
            ); ree;

        if (!ms2GroupMetadata->isValid()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK MS2 group metadata invalid"
                     << "group_id" << ms2GroupMetadata->groupId;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parseMetadataObject(
        const QJsonObject &jsonObject,
        const QString &sidecarRootPath,
        TimsbukIndexMetadata *metadata
        ) {

        ERR_INIT

        metadata->clear();

        e = parseRequiredStringField(
            jsonObject,
            JSON_FIELD_VERSION,
            QStringLiteral("metadata"),
            &metadata->version
            ); ree;
        e = parseRequiredStringField(
            jsonObject,
            JSON_FIELD_CREATED_AT,
            QStringLiteral("metadata"),
            &metadata->createdAt
            ); ree;

        const QJsonValue ms1PeaksValue = jsonObject.value(JSON_FIELD_MS1_PEAKS);
        if (!ms1PeaksValue.isObject()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK ms1_peaks must be an object";
            rrr(eFileError);
        }

        e = parsePeakGroupMetadata(
            ms1PeaksValue.toObject(),
            sidecarRootPath,
            JSON_FIELD_MS1_PEAKS,
            &metadata->ms1Peaks
            ); ree;

        const QJsonValue ms2WindowGroupsValue = jsonObject.value(JSON_FIELD_MS2_WINDOW_GROUPS);
        if (!ms2WindowGroupsValue.isArray()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK ms2_window_groups must be an array";
            rrr(eFileError);
        }

        const QJsonArray ms2WindowGroupsArray = ms2WindowGroupsValue.toArray();
        if (ms2WindowGroupsArray.isEmpty()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK ms2_window_groups is empty";
            rrr(eFileError);
        }

        QSet<TimsbukWindowGroupId> seenGroupIds;
        metadata->ms2WindowGroups.reserve(ms2WindowGroupsArray.size());
        for (const QJsonValue &groupValue : ms2WindowGroupsArray) {
            if (!groupValue.isObject()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK ms2_window_groups entry must be an object";
                rrr(eFileError);
            }

            TimsbukMs2GroupMetadata groupMetadata;
            e = parseMs2GroupMetadata(
                groupValue.toObject(),
                sidecarRootPath,
                &groupMetadata
                ); ree;

            if (seenGroupIds.contains(groupMetadata.groupId)) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK duplicate ms2_window_groups id"
                         << groupMetadata.groupId;
                rrr(eFileError);
            }

            seenGroupIds.insert(groupMetadata.groupId);
            metadata->ms2WindowGroups.push_back(groupMetadata);
        }

        if (!metadata->isValid()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK metadata invalid after parsing"
                     << "sidecar_root" << sidecarRootPath;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parseNumericRangeField(
        const QJsonObject &jsonObject,
        const QString &fieldName,
        const QString &context,
        float *lower,
        float *upper
        ) {

        ERR_INIT

        const QJsonValue fieldValue = jsonObject.value(fieldName);
        if (!fieldValue.isArray()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK isolation field must be an array"
                     << "context" << context
                     << "field" << fieldName;
            rrr(eFileError);
        }

        const QJsonArray jsonArray = fieldValue.toArray();
        constexpr int expectedSize = 2;
        if (jsonArray.size() != expectedSize || !jsonArray.at(0).isDouble() || !jsonArray.at(1).isDouble()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK isolation field must contain two numeric entries"
                     << "context" << context
                     << "field" << fieldName;
            rrr(eFileError);
        }

        *lower = static_cast<float>(jsonArray.at(0).toDouble());
        *upper = static_cast<float>(jsonArray.at(1).toDouble());
        if (*upper <= *lower) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK isolation field bounds invalid"
                     << "context" << context
                     << "field" << fieldName
                     << "lower" << *lower
                     << "upper" << *upper;
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err parseIsolationWindowSegments(
        const TimsbukMs2GroupMetadata &groupMetadata,
        QVector<TimsbukIsolationWindowSegment> *segments
        ) {

        ERR_INIT

        segments->clear();

        int windowIndex = 0;
        for (const QJsonValue &entryValue : groupMetadata.quadrupoleIsolation) {
            if (!entryValue.isObject()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK isolation entry must be an object"
                         << "group_id" << groupMetadata.groupId;
                rrr(eFileError);
            }

            const QJsonObject entryObject = entryValue.toObject();
            for (auto it = entryObject.constBegin(); it != entryObject.constEnd(); ++it) {
                if (!it.value().isObject()) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK isolation payload must be an object"
                             << "group_id" << groupMetadata.groupId
                             << "label" << it.key();
                    rrr(eFileError);
                }

                TimsbukIsolationWindowSegment segment;
                segment.groupId = groupMetadata.groupId;
                segment.windowIndex = windowIndex++;
                segment.windowLabel = it.key();
                segment.relativePath = groupMetadata.groupInfo.relativePath;

                const QString context = QStringLiteral("ms2_window_groups[%1].quadrupole_isolation[%2]")
                    .arg(groupMetadata.groupId)
                    .arg(segment.windowIndex);
                const QJsonObject payloadObject = it.value().toObject();

                e = parseNumericRangeField(
                    payloadObject,
                    JSON_FIELD_MZ,
                    context,
                    &segment.precursorMzLower,
                    &segment.precursorMzUpper
                    ); ree;
                e = parseNumericRangeField(
                    payloadObject,
                    JSON_FIELD_IM,
                    context,
                    &segment.ionMobilityLower,
                    &segment.ionMobilityUpper
                    ); ree;

                if (!segment.isValid()) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK isolation segment invalid after parsing"
                             << "context" << context;
                    rrr(eFileError);
                }

                segments->push_back(segment);
            }
        }

        e = ErrorUtils::isNotEmpty(*segments, eFileError); ree;

        ERR_RETURN
    }

    arrow::Status initArrowReaderBuilder(
        const QString &parquetFilePath,
        parquet::arrow::FileReaderBuilder *readerBuilder
        ) {

        arrow::MemoryPool *pool = arrow::default_memory_pool();

        auto readerProperties = parquet::ReaderProperties(pool);
        readerProperties.set_buffer_size(4096 * 4);
        readerProperties.enable_buffered_stream();

        auto arrowReaderProperties = parquet::ArrowReaderProperties();
        arrowReaderProperties.set_batch_size(128 * 1024);

        arrow::Status status = readerBuilder->OpenFile(
            parquetFilePath.toStdString(),
            /*memory_map=*/false,
            readerProperties
            );
        if (!status.ok()) {
            return status;
        }

        readerBuilder->memory_pool(pool);
        readerBuilder->properties(arrowReaderProperties);

        return status;
    }

    Err fieldIndexFromSchema(
        const std::shared_ptr<arrow::Schema> &schema,
        const std::string &fieldName,
        int *fieldIndex
        ) {

        ERR_INIT

        *fieldIndex = schema->GetFieldIndex(fieldName);
        if (*fieldIndex < 0) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK parquet column missing"
                     << QString::fromStdString(fieldName)
                     << QString::fromStdString(schema->ToString());
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err readFloatValue(
        const std::shared_ptr<arrow::Array> &array,
        int64_t rowIndex,
        float *value
        ) {

        ERR_INIT

        switch (array->type_id()) {
        case arrow::Type::FLOAT:
            *value = std::static_pointer_cast<arrow::FloatArray>(array)->Value(rowIndex);
            break;
        case arrow::Type::DOUBLE:
            *value = static_cast<float>(std::static_pointer_cast<arrow::DoubleArray>(array)->Value(rowIndex));
            break;
        default:
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK expected float-compatible parquet column"
                     << array->type()->ToString().c_str();
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err readIonMobilityValue(
        const std::shared_ptr<arrow::Array> &array,
        int64_t rowIndex,
        float *value
        ) {

        ERR_INIT

        switch (array->type_id()) {
        case arrow::Type::HALF_FLOAT: {
            const uint16_t bits = std::static_pointer_cast<arrow::HalfFloatArray>(array)->Value(rowIndex);
            *value = arrow::util::Float16::FromBits(bits).ToFloat();
            break;
        }
        case arrow::Type::FLOAT:
            *value = std::static_pointer_cast<arrow::FloatArray>(array)->Value(rowIndex);
            break;
        case arrow::Type::DOUBLE:
            *value = static_cast<float>(std::static_pointer_cast<arrow::DoubleArray>(array)->Value(rowIndex));
            break;
        default:
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK expected ion-mobility parquet column"
                     << array->type()->ToString().c_str();
            rrr(eFileError);
        }

        ERR_RETURN
    }

    Err readCycleIndexValue(
        const std::shared_ptr<arrow::Array> &array,
        int64_t rowIndex,
        TimsbukCycleIndex *cycleIndex
        ) {

        ERR_INIT

        switch (array->type_id()) {
        case arrow::Type::UINT32:
            *cycleIndex = std::static_pointer_cast<arrow::UInt32Array>(array)->Value(rowIndex);
            break;
        case arrow::Type::INT32: {
            const int32_t value = std::static_pointer_cast<arrow::Int32Array>(array)->Value(rowIndex);
            e = ErrorUtils::isTrue(value >= 0, eFileError); ree;
            *cycleIndex = static_cast<TimsbukCycleIndex>(value);
            break;
        }
        case arrow::Type::INT64: {
            const int64_t value = std::static_pointer_cast<arrow::Int64Array>(array)->Value(rowIndex);
            e = ErrorUtils::isTrue(0 <= value && value <= std::numeric_limits<quint32>::max(), eFileError); ree;
            *cycleIndex = static_cast<TimsbukCycleIndex>(value);
            break;
        }
        default:
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK expected cycle-index parquet column"
                     << array->type()->ToString().c_str();
            rrr(eFileError);
        }

        ERR_RETURN
    }

    int bestIsolationWindowIndex(
        const QVector<TimsbukIsolationWindowSegment> &segments,
        float ionMobility,
        bool *usedFallback = nullptr
        ) {

        if (usedFallback != nullptr) {
            *usedFallback = false;
        }

        for (int i = 0; i < segments.size(); ++i) {
            if (segments.at(i).containsIonMobility(ionMobility)) {
                return i;
            }
        }

        if (segments.isEmpty()) {
            return -1;
        }

        int bestIndex = 0;
        float bestDistance = std::abs(segments.first().ionMobilityCenter() - ionMobility);
        for (int i = 1; i < segments.size(); ++i) {
            const float distance = std::abs(segments.at(i).ionMobilityCenter() - ionMobility);
            if (distance < bestDistance) {
                bestDistance = distance;
                bestIndex = i;
            }
        }

        if (usedFallback != nullptr) {
            *usedFallback = true;
        }

        return bestIndex;
    }

    Err materializeMs2GroupScanStores(
        const QString &sidecarRootPath,
        TimsbukPendingMs2Group *pendingGroup
        ) {

        ERR_INIT

        e = parseIsolationWindowSegments(pendingGroup->metadata, &pendingGroup->isolationWindows); ree;
        e = ErrorUtils::isAboveThreshold(pendingGroup->cycleCount(), 1, ErrorUtilsParam::ExcludeThreshold, eFileError); ree;

        pendingGroup->scanStores.clear();
        pendingGroup->scanStores.resize(pendingGroup->isolationWindows.size() * pendingGroup->cycleCount());

        const QString parquetFilePath = cleanPath(
            QDir(sidecarRootPath).filePath(pendingGroup->metadata.groupInfo.relativePath)
            );

        parquet::arrow::FileReaderBuilder readerBuilder;
        arrow::Status status = initArrowReaderBuilder(parquetFilePath, &readerBuilder);
        if (!status.ok()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK parquet reader initialization failed"
                     << parquetFilePath
                     << QString::fromStdString(status.ToString());
            rrr(eFileError);
        }

        auto fileReaderResult = readerBuilder.Build();
        if (!fileReaderResult.ok()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK parquet reader build failed"
                     << parquetFilePath
                     << QString::fromStdString(fileReaderResult.status().ToString());
            rrr(eFileError);
        }

        std::unique_ptr<parquet::arrow::FileReader> fileReader = std::move(fileReaderResult).ValueOrDie();
        std::shared_ptr<arrow::RecordBatchReader> recordBatchReader;
        status = fileReader->GetRecordBatchReader(&recordBatchReader);
        if (!status.ok()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK parquet record-batch reader failed"
                     << parquetFilePath
                     << QString::fromStdString(status.ToString());
            rrr(eFileError);
        }

        int mzFieldIndex = -1;
        int intensityFieldIndex = -1;
        int ionMobilityFieldIndex = -1;
        int cycleIndexFieldIndex = -1;
        qint64 fallbackAssignedPeakCount = 0;
        qint64 totalPeakCount = 0;

        while (true) {
            std::shared_ptr<arrow::RecordBatch> recordBatch;
            status = recordBatchReader->ReadNext(&recordBatch);
            if (!status.ok()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMSBUK parquet batch read failed"
                         << parquetFilePath
                         << QString::fromStdString(status.ToString());
                rrr(eFileError);
            }

            if (!recordBatch) {
                break;
            }

            if (mzFieldIndex < 0) {
                const std::shared_ptr<arrow::Schema> schema = recordBatch->schema();
                e = fieldIndexFromSchema(schema, JSON_FIELD_MZ.toStdString(), &mzFieldIndex); ree;
                e = fieldIndexFromSchema(schema, QStringLiteral("intensity").toStdString(), &intensityFieldIndex); ree;
                e = fieldIndexFromSchema(schema, QStringLiteral("mobility_ook0").toStdString(), &ionMobilityFieldIndex); ree;
                e = fieldIndexFromSchema(schema, QStringLiteral("cycle_index").toStdString(), &cycleIndexFieldIndex); ree;
            }

            const std::shared_ptr<arrow::Array> mzArray = recordBatch->column(mzFieldIndex);
            const std::shared_ptr<arrow::Array> intensityArray = recordBatch->column(intensityFieldIndex);
            const std::shared_ptr<arrow::Array> ionMobilityArray = recordBatch->column(ionMobilityFieldIndex);
            const std::shared_ptr<arrow::Array> cycleIndexArray = recordBatch->column(cycleIndexFieldIndex);

            for (int64_t rowIndex = 0; rowIndex < recordBatch->num_rows(); ++rowIndex) {
                float mz = -1.0f;
                float intensity = -1.0f;
                float ionMobility = -1.0f;
                TimsbukCycleIndex cycleIndex = 0;

                e = readFloatValue(mzArray, rowIndex, &mz); ree;
                e = readFloatValue(intensityArray, rowIndex, &intensity); ree;
                e = readIonMobilityValue(ionMobilityArray, rowIndex, &ionMobility); ree;
                e = readCycleIndexValue(cycleIndexArray, rowIndex, &cycleIndex); ree;

                if (cycleIndex == 0 || cycleIndex >= static_cast<TimsbukCycleIndex>(pendingGroup->cycleCount())) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK cycle_index outside metadata range"
                             << "group_id" << pendingGroup->metadata.groupId
                             << "cycle_index" << cycleIndex
                             << "cycle_count" << pendingGroup->cycleCount();
                    rrr(eFileError);
                }

                bool usedFallback = false;
                const int windowIndex = bestIsolationWindowIndex(
                    pendingGroup->isolationWindows,
                    ionMobility,
                    &usedFallback
                    );
                if (windowIndex < 0) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMSBUK no isolation window available for peak"
                             << "group_id" << pendingGroup->metadata.groupId
                             << "ion_mobility" << ionMobility;
                    rrr(eFileError);
                }

                if (usedFallback) {
                    ++fallbackAssignedPeakCount;
                }

                pendingGroup->scanStore(windowIndex, cycleIndex)->scanPoints.push_back({mz, intensity});
                ++totalPeakCount;
            }
        }

        if (fallbackAssignedPeakCount > 0) {
            qWarning() << qPrintable(S_GLOBAL_TIMER.elapsed())
                       << "TIMSBUK peaks outside literal quadrupole isolation IM bounds were assigned to the nearest window"
                       << "group_id" << pendingGroup->metadata.groupId
                       << "fallback_assigned_peaks" << fallbackAssignedPeakCount
                       << "total_peaks" << totalPeakCount
                       << "relative_path" << pendingGroup->metadata.groupInfo.relativePath;
        }

        ERR_RETURN
    }

    Err appendPendingLogicalMs2Scans(
        const TimsbukPendingMs2Group &pendingGroup,
        QVector<TimsbukPendingLogicalMs2Scan> *pendingScans
        ) {

        ERR_INIT

        for (int windowIndex = 0; windowIndex < pendingGroup.isolationWindows.size(); ++windowIndex) {
            const TimsbukIsolationWindowSegment &window = pendingGroup.isolationWindows.at(windowIndex);
            for (TimsbukCycleIndex cycleIndex = 1;
                 cycleIndex < static_cast<TimsbukCycleIndex>(pendingGroup.cycleCount());
                 ++cycleIndex) {

                const TimsbukScanPointStore &pointStore = pendingGroup.scanStore(windowIndex, cycleIndex);
                if (pointStore.isEmpty()) {
                    continue;
                }

                TimsbukPendingLogicalMs2Scan pendingScan;
                pendingScan.isolationWindow = window;
                pendingScan.logicalScan.descriptor.msLevel = 2;
                pendingScan.logicalScan.descriptor.cycleIndex = cycleIndex;
                pendingScan.logicalScan.descriptor.windowGroupId = window.logicalWindowId();
                pendingScan.logicalScan.descriptor.scanTimeMilliseconds
                    = pendingGroup.metadata.groupInfo.cycleToRtMilliseconds.at(static_cast<int>(cycleIndex));
                pendingScan.logicalScan.descriptor.collisionEnergy = window.collisionEnergy;
                pendingScan.logicalScan.descriptor.precursorTargetMz = window.precursorTargetMz();
                pendingScan.logicalScan.descriptor.isoWindowLower = window.isolationWindowLower();
                pendingScan.logicalScan.descriptor.isoWindowUpper = window.isolationWindowUpper();
                pendingScan.logicalScan.descriptor.targetKey = window.targetKey();
                pendingScan.logicalScan.pointStore = pointStore;
                MsReaderBase::sortScanPoints(
                    ScanPointsSort::AscMz,
                    &pendingScan.logicalScan.pointStore.scanPoints
                    );

                pendingScans->push_back(pendingScan);
            }
        }

        ERR_RETURN
    }

    Err populateMs2ScansFromSidecar(
        const QString &sidecarRootPath,
        const TimsbukIndexMetadata &metadata,
        QMap<ScanNumber, MsScanInfo> *msScanInfo,
        QMap<ScanNumber, ScanPoints> *scanPoints,
        QMap<MzTargetKey, QVector<MsScanInfo*>> *mzTargetVsScanInfosPntrs,
        float *mzMs2Min,
        float *mzMs2Max
        ) {

        ERR_INIT

        msScanInfo->clear();
        scanPoints->clear();
        mzTargetVsScanInfosPntrs->clear();
        *mzMs2Min = std::numeric_limits<float>::max();
        *mzMs2Max = -1.0f;

        QVector<TimsbukPendingLogicalMs2Scan> pendingScans;
        int totalPossibleScanCount = 0;
        for (const TimsbukMs2GroupMetadata &groupMetadata : metadata.ms2WindowGroups) {
            QVector<TimsbukIsolationWindowSegment> isolationWindows;
            e = parseIsolationWindowSegments(groupMetadata, &isolationWindows); ree;
            totalPossibleScanCount += (groupMetadata.groupInfo.cycleToRtMilliseconds.size() - 1) * isolationWindows.size();
        }
        pendingScans.reserve(totalPossibleScanCount);

        for (const TimsbukMs2GroupMetadata &groupMetadata : metadata.ms2WindowGroups) {
            TimsbukPendingMs2Group pendingGroup;
            pendingGroup.metadata = groupMetadata;

            e = materializeMs2GroupScanStores(sidecarRootPath, &pendingGroup); ree;
            e = appendPendingLogicalMs2Scans(pendingGroup, &pendingScans); ree;
        }

        e = ErrorUtils::isNotEmpty(pendingScans, eFileError); ree;

        std::stable_sort(
            pendingScans.begin(),
            pendingScans.end(),
            [](const TimsbukPendingLogicalMs2Scan &left, const TimsbukPendingLogicalMs2Scan &right) {
                if (!MathUtils::tSame(
                        left.logicalScan.descriptor.scanTimeMilliseconds,
                        right.logicalScan.descriptor.scanTimeMilliseconds
                        )) {
                    return left.logicalScan.descriptor.scanTimeMilliseconds
                         < right.logicalScan.descriptor.scanTimeMilliseconds;
                }

                if (left.logicalScan.descriptor.windowGroupId != right.logicalScan.descriptor.windowGroupId) {
                    return left.logicalScan.descriptor.windowGroupId
                         < right.logicalScan.descriptor.windowGroupId;
                }

                return left.logicalScan.descriptor.precursorTargetMz
                     < right.logicalScan.descriptor.precursorTargetMz;
            }
            );

        ScanNumber scanNumber = 1;
        qint64 totalPointCount = 0;
        for (const TimsbukPendingLogicalMs2Scan &pendingScan : pendingScans) {
            const ScanPoints &logicalScanPoints = pendingScan.logicalScan.pointStore.scanPoints;
            if (logicalScanPoints.isEmpty()) {
                continue;
            }

            MsScanInfo logicalMsScanInfo;
            logicalMsScanInfo.msLevel = 2;
            logicalMsScanInfo.scanNumber = scanNumber;
            logicalMsScanInfo.scanTime = pendingScan.logicalScan.descriptor.scanTimeMilliseconds
                * TIMSBUK_SCAN_TIME_MILLISECONDS_TO_MINUTES;
            logicalMsScanInfo.collisionEnergy = pendingScan.logicalScan.descriptor.collisionEnergy;
            logicalMsScanInfo.precursorTargetMz = pendingScan.logicalScan.descriptor.precursorTargetMz;
            logicalMsScanInfo.isoWindowLower = pendingScan.logicalScan.descriptor.isoWindowLower;
            logicalMsScanInfo.isoWindowUpper = pendingScan.logicalScan.descriptor.isoWindowUpper;
            logicalMsScanInfo.ionMobilityDriftTime = pendingScan.isolationWindow.ionMobilityCenter();
            logicalMsScanInfo.ionMobilityWindowLower = pendingScan.isolationWindow.ionMobilityLower;
            logicalMsScanInfo.ionMobilityWindowUpper = pendingScan.isolationWindow.ionMobilityUpper;
            logicalMsScanInfo.nativeFrameNumber = static_cast<int>(pendingScan.logicalScan.descriptor.cycleIndex);
            logicalMsScanInfo.nativeScanNumber = pendingScan.logicalScan.descriptor.windowGroupId;

            msScanInfo->insert(scanNumber, logicalMsScanInfo);
            scanPoints->insert(scanNumber, logicalScanPoints);

            *mzMs2Min = std::min(*mzMs2Min, logicalScanPoints.front().x());
            *mzMs2Max = std::max(*mzMs2Max, logicalScanPoints.back().x());
            totalPointCount += logicalScanPoints.size();
            ++scanNumber;
        }

        e = ErrorUtils::isNotEmpty(*msScanInfo, eFileError); ree;
        e = ErrorUtils::isNotEmpty(*scanPoints, eFileError); ree;

        for (auto it = msScanInfo->begin(); it != msScanInfo->end(); ++it) {
            (*mzTargetVsScanInfosPntrs)[it.value().targetKey()].push_back(&it.value());
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "MsReaderTimsbukIndex materialized MS2 logical scans"
                 << "scan_count" << msScanInfo->size()
                 << "target_count" << mzTargetVsScanInfosPntrs->size()
                 << "point_count" << totalPointCount;

        ERR_RETURN
    }

}//namespace

bool MsReaderTimsbukIndex::isDirectIndexRootPath(const QString &filePath) {
    if (filePath.isEmpty()) {
        return false;
    }

    const QString normalizedPath = cleanPath(filePath);
    const QFileInfo fileInfo(normalizedPath);
    return fileInfo.exists()
        && fileInfo.isDir()
        && directoryContainsMetadataFile(normalizedPath);
}

Err MsReaderTimsbukIndex::resolveInputPath(
    const QString &inputPath,
    QString *sidecarRootPath,
    QString *sourceBrukerDirectoryPath
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(inputPath); ree;

    const QString normalizedInputPath = cleanPath(inputPath);
    const QFileInfo inputInfo(normalizedInputPath);
    if (!inputInfo.exists()) {
        if (normalizedInputPath.endsWith(QStringLiteral(".idx"), Qt::CaseInsensitive)) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMSBUK sidecar path missing"
                     << normalizedInputPath;
            rrr(eFileError);
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "Unsupported TIMSBUK input path"
                 << normalizedInputPath
                 << "Expected Bruker .d directory or sidecar root";
        rrr(eFileIncorrectTypeError);
    }

    if (!inputInfo.isDir()) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "Unsupported TIMSBUK input path type"
                 << normalizedInputPath;
        rrr(eFileIncorrectTypeError);
    }

    if (hasBrukerDirectorySuffix(inputInfo)) {
        const QString normalizedSidecarRootPath = normalizedInputPath + QStringLiteral(".idx");
        e = validateSidecarRootPath(normalizedSidecarRootPath, normalizedInputPath); ree;

        if (sidecarRootPath != nullptr) {
            *sidecarRootPath = normalizedSidecarRootPath;
        }
        if (sourceBrukerDirectoryPath != nullptr) {
            *sourceBrukerDirectoryPath = normalizedInputPath;
        }

        ERR_RETURN
    }

    if (isDirectIndexRootPath(normalizedInputPath)) {
        if (sidecarRootPath != nullptr) {
            *sidecarRootPath = normalizedInputPath;
        }
        if (sourceBrukerDirectoryPath != nullptr) {
            *sourceBrukerDirectoryPath = inferSourceBrukerDirectoryPath(normalizedInputPath);
        }

        ERR_RETURN
    }

    if (hasIndexDirectorySuffix(inputInfo)) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMSBUK sidecar root missing metadata.json"
                 << normalizedInputPath;
        rrr(eFileError);
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Unsupported TIMSBUK input path"
             << normalizedInputPath
             << "Expected Bruker .d directory or sidecar root";
    rrr(eFileIncorrectTypeError);
}

Err MsReaderTimsbukIndex::loadMetadata(
    const QString &sidecarRootPath,
    TimsbukIndexMetadata *metadata
    ) {

    ERR_INIT

    const QString metadataFilePath = metadataFilePathForRoot(sidecarRootPath);

    QFile metadataFile(metadataFilePath);
    if (!metadataFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMSBUK metadata.json could not be opened"
                 << metadataFilePath;
        rrr(eFileError);
    }

    const QByteArray metadataBytes = metadataFile.readAll();
    metadataFile.close();
    e = ErrorUtils::isNotEmpty(metadataBytes, eFileError); ree;

    QJsonParseError parseError;
    const QJsonDocument jsonDocument = QJsonDocument::fromJson(metadataBytes, &parseError);
    if (parseError.error != QJsonParseError::NoError || !jsonDocument.isObject()) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMSBUK metadata.json parse failed"
                 << metadataFilePath
                 << parseError.errorString();
        rrr(eFileError);
    }

    e = parseMetadataObject(
        jsonDocument.object(),
        sidecarRootPath,
        metadata
        ); ree;

    ERR_RETURN
}

Err MsReaderTimsbukIndex::openFile(const QString &filePath) {

    ERR_INIT

    QString sidecarRootPath;
    QString sourceBrukerDirectoryPath;
    TimsbukIndexMetadata metadata;

    e = closeFile(); ree;
    e = resolveInputPath(
        filePath,
        &sidecarRootPath,
        &sourceBrukerDirectoryPath
        ); ree;
    e = loadMetadata(sidecarRootPath, &metadata); ree;

    m_sidecarRootPath = sidecarRootPath;
    m_sourceBrukerDirectoryPath = sourceBrukerDirectoryPath;
    m_metadata = metadata;
    m_metadataFilePath = metadataFilePathForRoot(m_sidecarRootPath);
    m_filePath = m_sidecarRootPath;
    m_isTIMS = false;

    e = populateMs2ScansFromSidecar(
        m_sidecarRootPath,
        m_metadata,
        &m_msScanInfo,
        &m_scanPoints,
        &m_mzTargetVsScanInfosPntrs,
        &m_mzMs2Min,
        &m_mzMs2Max
        ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex sidecar reader path active"
             << "input_path" << filePath
             << "sidecar_root" << m_sidecarRootPath
             << "source_bruker" << m_sourceBrukerDirectoryPath
             << "metadata_version" << m_metadata.version
             << "ms2_group_count" << m_metadata.ms2WindowGroups.size()
             << "scan_count" << m_msScanInfo.size();
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex metadata loaded and MS2 scans materialized";

    ERR_RETURN
}

Err MsReaderTimsbukIndex::openFile(
    const QString &filePath,
    const QString &columnToFilterBy,
    const QPair<double, double> &filterRange
    ) {

    Q_UNUSED(columnToFilterBy)
    Q_UNUSED(filterRange)

    return openFile(filePath);
}

Err MsReaderTimsbukIndex::getMzTargetScanPoints(
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

Err MsReaderTimsbukIndex::closeFile() {
    ERR_INIT

    m_sidecarRootPath.clear();
    m_sourceBrukerDirectoryPath.clear();
    m_metadataFilePath.clear();
    m_metadata.clear();
    m_mzTargetVsScanInfosPntrs.clear();
    m_frameIndexVsDriftTime.clear();
    m_filePath.clear();
    m_isTIMS = false;
    m_mzMs1Min = std::numeric_limits<float>::max();
    m_mzMs1Max = -1.0f;
    m_mzMs2Min = std::numeric_limits<float>::max();
    m_mzMs2Max = -1.0f;

    e = MsReaderBase::closeFile(); ree;

    ERR_RETURN
}
