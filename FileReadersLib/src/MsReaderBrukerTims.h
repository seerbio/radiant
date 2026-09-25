//
// Created by andrewnichols on 10/8/24.
//

#ifndef MSREADERBRUKERTIMS_H
#define MSREADERBRUKERTIMS_H

#include "FileReadersLib_Exports.h"

#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "SqlUtils.h"

using namespace Error;

struct TimsMS2WindowsInfo {
    int windowGroup = -1;
    int scanNumberBegin = -1;
    int scanNumberEnd = -1;
    float isolationMz = -1.0;
    float isolationWidth = -1.0;
    float collisionEnergy = -1.0;
};

struct TimsFrameInfo {
    int frameId = -1;
    float scanTime = -1.0;
    int msmsType = -1;
    int numScans = -1;
    int windowGroup = -1;
};

class FILEREADERSLIB_EXPORTS MsReaderBrukerTims : public MsReaderBase {

public:

    MsReaderBrukerTims() = default;

    ~MsReaderBrukerTims() override = default;

    Err openFile(const QString &filePath) override;

    Err openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) override;

    Err closeFile() override;

    Err writeFrame(
        const QString &filePath,
        float scanTime,
        int msLevel
        );

private:

    Err openFile(
        const QString &filePath,
        bool filterByScanTime,
        const QPair<double, double> &scanTimeFilterRange
        );

    QHash<MzTargetKey, TimsMS2WindowsInfo> m_mzTargetVsTimsMs2WindowsInfos;
    QHash<int, QVector<TimsMS2WindowsInfo>> m_windowGroupIndexVsTimsMs2WindowsInfoses;
    std::vector<TimsFrameInfo> m_timsFramesInfos;
};



#endif //MSREADERBRUKERTIMS_H
