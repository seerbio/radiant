//
// Created by Codex on 9/8/26.
//

#ifndef TIMSREADERRUNCATALOG_H
#define TIMSREADERRUNCATALOG_H

#include "Error.h"
#include "GlobalSettings.h"

#include <QMap>
#include <QString>
#include <QVector>

#include <cstdint>

using namespace Error;

struct TimsreaderIsolationWindowInfo {
    uint32_t isolationWindowId = 0;
    int32_t windowGroupId = -1;
    uint16_t scanStartIndex = 0;
    uint16_t scanEndIndexExclusive = 0;
    float isolationTargetMz = -1.0f;
    float isolationLowerOffset = -1.0f;
    float isolationUpperOffset = -1.0f;

    [[nodiscard]] bool isValid() const;
    [[nodiscard]] MzTargetKey targetKey() const;
};

class TimsreaderRunCatalog {
public:
    TimsreaderRunCatalog() = default;

    Err open(const QString &brukerPath);
    void close();

    [[nodiscard]] bool isInit() const;
    [[nodiscard]] QString filePath() const;
    [[nodiscard]] QVector<TimsreaderIsolationWindowInfo> isolationWindows() const;
    [[nodiscard]] QVector<uint32_t> isolationWindowIdsForTargetKey(
        const MzTargetKey &targetKey
        ) const;

private:
    QString m_filePath;
    QVector<TimsreaderIsolationWindowInfo> m_isolationWindows;
    QMap<MzTargetKey, QVector<uint32_t>> m_targetKeyToIsolationWindowIds;
};

#endif // TIMSREADERRUNCATALOG_H
