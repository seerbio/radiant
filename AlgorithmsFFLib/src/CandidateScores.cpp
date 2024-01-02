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
