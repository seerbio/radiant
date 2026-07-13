//
// Created by Codex on 7/13/26.
//

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

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex sidecar reader path active"
             << "input_path" << filePath
             << "sidecar_root" << m_sidecarRootPath
             << "source_bruker" << m_sourceBrukerDirectoryPath
             << "metadata_version" << m_metadata.version
             << "ms2_group_count" << m_metadata.ms2WindowGroups.size();
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex metadata loaded; scan materialization pending";

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
