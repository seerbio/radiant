//
// Created by Codex on 9/8/26.
//

#include "TimsreaderRunCatalog.h"

#include "ErrorUtils.h"
#include "MsReaderBase.h"

#include <QDebug>

#ifdef PYTHIA_HAVE_TIMSREADER
#include "timsreader.hpp"
#endif

namespace {

#ifdef PYTHIA_HAVE_TIMSREADER
TimsreaderIsolationWindowInfo makeIsolationWindowInfo(
    const timsreader::IsolationWindowInfo &window
    ) {

    TimsreaderIsolationWindowInfo info;
    info.isolationWindowId = window.isolation_window_id;
    info.windowGroupId = window.window_group_id;
    info.scanStartIndex = window.scan_start_index;
    info.scanEndIndexExclusive = window.scan_end_index_exclusive;
    info.isolationTargetMz = window.isolation_target_mz;
    info.isolationLowerOffset = window.isolation_lower_offset;
    info.isolationUpperOffset = window.isolation_upper_offset;
    return info;
}
#endif

}

bool TimsreaderIsolationWindowInfo::isValid() const {
    return isolationTargetMz > 0.0f
        && isolationLowerOffset >= 0.0f
        && isolationUpperOffset >= 0.0f;
}

MzTargetKey TimsreaderIsolationWindowInfo::targetKey() const {
    return MsScanInfo::targetKey(
        isolationTargetMz - isolationLowerOffset,
        isolationTargetMz + isolationUpperOffset
        );
}

Err TimsreaderRunCatalog::open(const QString &brukerPath) {

    ERR_INIT

    close();
    e = ErrorUtils::fileExists(brukerPath); ree;

#ifdef PYTHIA_HAVE_TIMSREADER
    try {
        timsreader::Run run(brukerPath.toStdString());
        const std::vector<timsreader::IsolationWindowInfo> windows = run.isolation_windows();

        for (const timsreader::IsolationWindowInfo &window : windows) {
            const TimsreaderIsolationWindowInfo info = makeIsolationWindowInfo(window);
            e = ErrorUtils::isTrue(info.isValid(), eFileError); ree;
            m_isolationWindows.push_back(info);
            m_targetKeyToIsolationWindowIds[info.targetKey()].push_back(info.isolationWindowId);
        }
    }
    catch (const std::exception &ex) {
        qDebug() << "timsreader failed to open Bruker run" << brukerPath << ex.what();
        close();
        rrr(eFileError);
    }
#else
    Q_UNUSED(brukerPath)
    rrr(eFunctionNotImplemented);
#endif

    m_filePath = brukerPath;

    ERR_RETURN
}

void TimsreaderRunCatalog::close() {
    m_filePath.clear();
    m_isolationWindows.clear();
    m_targetKeyToIsolationWindowIds.clear();
}

bool TimsreaderRunCatalog::isInit() const {
    return !m_filePath.isEmpty();
}

QString TimsreaderRunCatalog::filePath() const {
    return m_filePath;
}

QVector<TimsreaderIsolationWindowInfo> TimsreaderRunCatalog::isolationWindows() const {
    return m_isolationWindows;
}

QVector<uint32_t> TimsreaderRunCatalog::isolationWindowIdsForTargetKey(
    const MzTargetKey &targetKey
    ) const {
    return m_targetKeyToIsolationWindowIds.value(targetKey);
}
