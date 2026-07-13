//
// Created by Codex on 7/13/26.
//

#include "MsReaderTimsbukIndex.h"

#include "ErrorUtils.h"
#include "StringUtils.h"

#include <QDebug>
#include <QDir>
#include <QFileInfo>


namespace {

    const QString TIMSBUK_INDEX_SUFFIX = QStringLiteral("idx");
    const QString TIMSBUK_METADATA_FILE_NAME = QStringLiteral("metadata.json");
    const QString TIMSBUK_BRUKER_INDEX_SUFFIX = QStringLiteral(".d.idx");

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

Err MsReaderTimsbukIndex::openFile(const QString &filePath) {

    ERR_INIT

    e = closeFile(); ree;
    e = resolveInputPath(
        filePath,
        &m_sidecarRootPath,
        &m_sourceBrukerDirectoryPath
        ); ree;

    m_filePath = m_sidecarRootPath;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex sidecar reader path active"
             << "input_path" << filePath
             << "sidecar_root" << m_sidecarRootPath
             << "source_bruker" << m_sourceBrukerDirectoryPath;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex parsing not implemented yet";

    rrr(eFunctionNotImplemented);
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
    m_sidecarRootPath.clear();
    m_sourceBrukerDirectoryPath.clear();
    return MsReaderBase::closeFile();
}
