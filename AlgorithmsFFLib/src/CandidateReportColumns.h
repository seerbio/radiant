#ifndef RADIANT_CANDIDATE_REPORT_COLUMNS_H
#define RADIANT_CANDIDATE_REPORT_COLUMNS_H

#include "CandidateBundleIO.h"

namespace CandidateReportColumns {
// Uses the same field definitions, types and Parquet settings as the row writer.
ALGORITHMSFFLIB_EXPORTS Error::Err write(
    const QVector<CandidateScoresReaderRow> &rows, const QString &path);

// written=false requests the legacy conversion path for noncanonical schemas.
// Selection and probability ranges have already been validated by the caller.
ALGORITHMSFFLIB_EXPORTS Error::Err tryWriteCombined(
    const QVector<CandidateBundleIO::View> &views,
    const CandidatePoolRescorer::Result &result,
    const QVector<double> &probabilities, const QVector<int> &destinations,
    const QString &path, int familyFolds, bool *written);
}

#endif
