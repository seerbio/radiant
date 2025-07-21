//
// Created by anichols on 10/20/23.
//

#include "CandidateScores.h"

#include "BiophysicalCalcs.h"


CandidateScores::CandidateScores(TargetDecoyCandidatePair *tdcp, bool isDecoy) {
    targetDecoyCandidatePair = tdcp;
    m_isDecoy = isDecoy;
}

void CandidateScores::initFeaturesArray() {
    featuresArray = QVector<float>(FeaturesSize, 0.0f);
    featuresArray.reserve(FeaturesSize);
}

bool CandidateScores::isDecoy() const {
    return m_isDecoy;
}

PeptideString CandidateScores::peptideString() const {
    return peptideStringWithMods().removeUniModChars();
}

PeptideStringWithMods CandidateScores::peptideStringWithMods() const {
    if (isDecoy()) {
        // Return the decoy sequence
        if (targetDecoyCandidatePair->isDecoy()) {
            // If the decoy is a library decoy, return the sequence from the library
            return targetDecoyCandidatePair->peptideStringWithMods();
        } else {
            // Otherwise, mutate the target sequence to get the decoy sequence
            return AminoAcids::mutatePenultimatePeptideResidues(targetDecoyCandidatePair->peptideStringWithMods());
        }
    } else {
        // Return the target sequence
        if (targetDecoyCandidatePair->isDecoy()) {
            // If the decoy is a library decoy, mutate the target sequence to get the decoy sequence
            return AminoAcids::mutatePenultimatePeptideResidues(targetDecoyCandidatePair->peptideStringWithMods());
        } else {
            // Simply return the target sequence
            return targetDecoyCandidatePair->peptideStringWithMods();
        }
    }
}

PeptideStringWithMods CandidateScores::peptideStringWithModsDecoyOrigin() const {
    if (targetDecoyCandidatePair->isDecoy()) {
        // If the decoy is a library decoy, return a mutated sequence; this is probably unimportant, but it does
        // mean that in most cases library decoys can not be correctly reannotated using this column!!!
        return AminoAcids::mutatePenultimatePeptideResidues(targetDecoyCandidatePair->peptideStringWithMods());
    }

    return targetDecoyCandidatePair->peptideStringWithMods();
}

float CandidateScores::mz() const {
    return targetDecoyCandidatePair->mz(isDecoy());
}

int CandidateScores::charge() const {
    return targetDecoyCandidatePair->charge();
}

float CandidateScores::mass() const {
    return targetDecoyCandidatePair->mass();
}

float CandidateScores::iRt() const {
    return targetDecoyCandidatePair->iRt();
}

float CandidateScores::iIM() const {
    return targetDecoyCandidatePair->iIM();
}

QVector<MS2Ion> CandidateScores::ms2Ions() const {
    if (isDecoy()) {
        // Return the ions used for scoring the decoy of this pair
        if (targetDecoyCandidatePair->isDecoy()) {
            // Use ions from the library
            return targetDecoyCandidatePair->ms2IonsTarget();
        } else {
            // Use shuffled ions
            return targetDecoyCandidatePair->ms2IonsDecoy();
        }
    } else {
        // Return the ions used for scoring the target of this pair
        if (targetDecoyCandidatePair->isDecoy()) {
            // Use shuffled ions
            return targetDecoyCandidatePair->ms2IonsDecoy();
        } else {
            // Use ions from the library
            return targetDecoyCandidatePair->ms2IonsTarget();
        }
    }
}

int CandidateScores::totalFragmentCount() const {
    return targetDecoyCandidatePair->totalFragmentCount();
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

void CandidateScores::printFeatures(const QVector<Features> &featuresToPrint) {

    QVector<float> vec;
    for (const Features &f : featuresToPrint) {
        vec.push_back(featuresArray[f]);
    }
    qDebug() << vec;
}
