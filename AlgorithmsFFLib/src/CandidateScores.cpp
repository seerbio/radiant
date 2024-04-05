//
// Created by anichols on 10/20/23.
//

#include "CandidateScores.h"

#include "BiophysicalCalcs.h"


void CandidateScores::initFeaturesArray() {
    featuresArray = QVector<float>(Features::FeaturesSize, 0.0f);
    featuresArray.reserve(Features::FeaturesSize);
}

QVector<float>* CandidateScores::featuresArrayRef() {
    return &featuresArray;
}

QVector<float> CandidateScores::selectFeaturesArrayFeatures(const QVector<Features> &enumFeatures) {

    QVector<float> selectedFeatures(enumFeatures.size());

    for (int i = 0; i < enumFeatures.size(); i++) {
        selectedFeatures[i] = featuresArray[enumFeatures.at(i)];
    }

    return selectedFeatures;
}

QVector<float> CandidateScores::selectFeaturesArrayFeatures(
        const QVector<float> &featureArray,
        const QVector<Features> &enumFeatures
        ) {
    QVector<float> selectedFeatures(enumFeatures.size());

    for (int i = 0; i < enumFeatures.size(); i++) {
        selectedFeatures[i] = featureArray[enumFeatures.at(i)];
    }

    return selectedFeatures;

}
