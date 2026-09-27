#ifndef RADIANT_CANDIDATE_POOL_RESCORER_H
#define RADIANT_CANDIDATE_POOL_RESCORER_H

#include "CandidatePoolSelection.h"
#include "PeptideFamilyNeuralNet.h"

// A complete native scoring operation over captured candidate views. The
// caller supplies aligned raw rows and owns file loading and report writing.
class ALGORITHMSFFLIB_EXPORTS CandidatePoolRescorer {
public:
    struct Run {
        QString name;
        QVector<CandidatePoolSelection::Identity> identities;
        QVector<FragmentCompetition::Evidence> evidence;
        QVector<QVector<float>> rawFeatures;
    };
    struct Settings {
        int originalCandidateLimit = 100000;
        int expandedCandidateLimit = 200000;
        int minimumSharedFragments = 2;
        PeptideFamilyNeuralNet::Settings neuralNet;
    };
    struct Result {
        QString name;
        // Indices refer to the caller's original run arrays.
        QVector<int> selectedInputIndices;
        QVector<float> originalPerRunProbability;
        QVector<float> originalPooledProbability;
        QVector<float> expandedPooledProbability;
        QVector<double> combinedProbability;
        // Confidence indices also refer to the original caller arrays.
        CandidatePoolSelection::Confidence confidence;
    };

    static double equalLogit(double first, double second);
    static Error::Err score(const QVector<Run> &runs, const Settings &settings,
                            QVector<Result> *results);
};

#endif
