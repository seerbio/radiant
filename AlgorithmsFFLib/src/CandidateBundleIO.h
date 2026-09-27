#ifndef RADIANT_CANDIDATE_BUNDLE_IO_H
#define RADIANT_CANDIDATE_BUNDLE_IO_H

#include "CandidatePoolRescorer.h"
#include <QJsonObject>

class ALGORITHMSFFLIB_EXPORTS CandidateBundleIO {
public:
    struct View {
        CandidatePoolRescorer::Run candidates;
        QJsonObject provenance;
        QString directory;
        QString reportSha256;
        QString manifestSha256;
        QJsonObject fileSha256;
    };

    // Rows must be in the native within-view LDA order, before NN transforms.
    // Writes a new directory and creates the completion manifest last.
    static Error::Err write(const QVector<CandidateScores*> &rows, const QString &runName,
                            const QJsonObject &provenance, const QString &directory);
    // File hashes, the exact feature schema and all row counts are checked.
    // Protein annotations stay in the separate report; the scorer does not
    // load that report. Output is replaced only on success.
    static Error::Err read(const QString &directory, View *view);
    // Results refer to concatenated views in the supplied order. Keeps all
    // original evidence columns and annotations, replacing precursor scores
    // and confidence. Peptide/protein confidence is reset to unassigned.
    static Error::Err writeReport(const QVector<View> &views,
                                 const CandidatePoolRescorer::Result &result,
                                 const QString &path, int familyFolds = 3);
    static QString fileHash(const QString &path);
};

#endif
