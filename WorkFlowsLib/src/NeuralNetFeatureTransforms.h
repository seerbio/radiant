#ifndef RADIANT_NEURAL_NET_FEATURE_TRANSFORMS_H
#define RADIANT_NEURAL_NET_FEATURE_TRANSFORMS_H

#include "CandidateScores.h"

#include <algorithm>
#include <cmath>

namespace NeuralNetFeatureTransforms {

inline void logIntensityFeatures(
    const QVector<Features> &features,
    QVector<float> *values
    ) {

    Q_ASSERT(values != nullptr);
    Q_ASSERT(features.size() == values->size());

    for (int index = 0; index < features.size(); ++index) {
        switch (features.at(index)) {
        case Ms1IntensityFound100:
        case Ms1IntensityFound45:
        case Ms1IntensityFoundPreMono:
        case Ms1IntensityFoundIso1:
        case Ms1IntensityFoundIso2:
        case Ms1IntensityFoundApex100IM:
            (*values)[index] = std::log1p(std::max(0.0f, (*values)[index]));
            break;
        default:
            break;
        }
    }
}

} // namespace NeuralNetFeatureTransforms

#endif // RADIANT_NEURAL_NET_FEATURE_TRANSFORMS_H
