#ifndef RADIANT_NEURAL_NET_FEATURE_TRANSFORMS_H
#define RADIANT_NEURAL_NET_FEATURE_TRANSFORMS_H

#include "CandidateScores.h"

#include <QCryptographicHash>
#include <QRegularExpression>
#include <QString>
#include <QVector>

#include <algorithm>
#include <cmath>

namespace NeuralNetFeatureTransforms {

inline QString peptideFamily(QString originSequence) {
    originSequence.remove(QRegularExpression("\\([^)]*\\)"));
    originSequence.remove('_');
    originSequence.replace('I', 'L');
    return originSequence;
}

// Use the target origin for both members of a target/decoy pair. The stable
// hash keeps charges, modifications and I/L equivalents in the same fold.
inline quint32 peptideFamilyHash(const QString &originSequence) {
    return QCryptographicHash::hash(peptideFamily(originSequence).toUtf8(),
                                   QCryptographicHash::Sha256).toHex().left(8).toUInt(nullptr, 16);
}

inline void logIntensities(const QVector<Features> &features, QVector<float> *values) {
    Q_ASSERT(features.size() == values->size());
    for (int index = 0; index < features.size(); ++index) {
        const auto feature = features[index];
        if (feature == TotalIntensityRaw || feature == Ms1IntensityFound100
            || feature == Ms1IntensityFound45 || feature == Ms1IntensityFoundPreMono
            || feature == Ms1IntensityFoundIso1 || feature == Ms1IntensityFoundIso2
            || feature == Ms1IntensityFoundApex100IM) {
            (*values)[index] = std::log1p(std::max(0.0f, (*values)[index]));
        }
    }
}

} // namespace NeuralNetFeatureTransforms

#endif
