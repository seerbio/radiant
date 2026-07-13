//
// Created by Codex on 7/13/26.
//

#ifndef MSREADERTIMSBUKINDEX_H
#define MSREADERTIMSBUKINDEX_H

#include "FileReadersLib_Exports.h"

#include "MsReaderBase.h"
#include "TimsbukIndexTypes.h"

using namespace Error;

class FILEREADERSLIB_EXPORTS MsReaderTimsbukIndex : public MsReaderBase {

public:

    MsReaderTimsbukIndex() = default;
    ~MsReaderTimsbukIndex() override = default;

    [[nodiscard]] static bool isDirectIndexRootPath(const QString &filePath);

    static Err resolveInputPath(
        const QString &inputPath,
        QString *sidecarRootPath,
        QString *sourceBrukerDirectoryPath
        );

    static Err loadMetadata(
        const QString &sidecarRootPath,
        TimsbukIndexMetadata *metadata
        );

    Err openFile(const QString &filePath) override;

    Err openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) override;

    Err getMzTargetScanPoints(
        const MzTargetKey &targetKey,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
        ) override;

    Err closeFile() override;

private:

    QString m_sidecarRootPath;
    QString m_sourceBrukerDirectoryPath;
    QString m_metadataFilePath;
    TimsbukIndexMetadata m_metadata;
};

#endif // MSREADERTIMSBUKINDEX_H
