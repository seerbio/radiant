//
// Created by Codex on 9/8/26.
//

#ifndef MSREADERTIMSREADER_H
#define MSREADERTIMSREADER_H

#include "FileReadersLib_Exports.h"

#include "ImHandlingMode.h"
#include "MsReaderBase.h"
#include "TimsbukIndexTypes.h"

#include <QScopedPointer>

class FILEREADERSLIB_EXPORTS MsReaderTimsreader : public MsReaderBase {
public:
    explicit MsReaderTimsreader(
        ImHandlingMode imHandlingMode = ImHandlingMode::Centroid,
        int threadCount = 1
        );
    ~MsReaderTimsreader() override;

    Err openFile(const QString &filePath) override;

    Err openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) override;

    Err restrictScanTimeRange(ScanTime scanTimeMin, ScanTime scanTimeMax) override;

    Err getMzTargetScanPoints(
        const MzTargetKey &targetKey,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
        ) override;

    Err getMzTargetAlignedPointData(
        const MzTargetKey &targetKey,
        QMap<ScanNumber, const TimsbukAlignedPointData*> *scanNumberVsAlignedPointData
        ) const override;

    const TimsbukAlignedPointData *alignedPointDataPntr(ScanNumber scanNumber) const override;

    Err closeFile() override;

private:
    Q_DISABLE_COPY(MsReaderTimsreader)

    class Private;
    QScopedPointer<Private> d_ptr;
};

#endif // MSREADERTIMSREADER_H
