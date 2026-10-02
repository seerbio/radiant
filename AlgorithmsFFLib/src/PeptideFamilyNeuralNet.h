#ifndef RADIANT_PEPTIDE_FAMILY_NEURAL_NET_H
#define RADIANT_PEPTIDE_FAMILY_NEURAL_NET_H

#include "CandidateScores.h"
#include "Error.h"

#include <QString>
#include <QVector>

// Stateless scoring of a caller-selected candidate pool. All runs and views
// of an origin peptide must use the same family hash. Identity annotations,
// species, proteins and reference identifications are deliberately absent.
class ALGORITHMSFFLIB_EXPORTS PeptideFamilyNeuralNet {
public:
    struct Settings {
        int folds = 3;
        int networks = 4;
        int epochs = 24;
        int threads = 3;
        int seed = 666;
        double learningRate = .003;
        double nodesFraction = .5;
        float focalLossGamma = 0.0f;
        bool logIntensities = true;
        bool shuffleEachEpoch = true;
        // Empty means one scale fitted across the complete selected pool.
        QVector<int> normalizationGroups;
    };

    struct Predictions {
        QVector<float> perNetwork; // Row-major [candidate][network].
        QVector<float> meanDecoyProbability;
        QVector<quint32> heldOutFold;
    };

    static quint32 familyHash(QString originSequence);
    static QVector<Features> nonTimsFeatures();

    // Output is replaced only on success. The caller owns candidate ranking,
    // view selection, physical competition and final confidence calculation.
    static Error::Err score(
        const QVector<QVector<float>> &rawFeatures,
        const QVector<quint32> &decoyLabels,
        const QVector<quint32> &familyHashes,
        const QVector<Features> &features,
        const Settings &settings,
        Predictions *predictions);
};

#endif
