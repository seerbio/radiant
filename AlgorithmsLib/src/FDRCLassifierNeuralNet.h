//
// Created by anichols on 9/26/23.
//

#ifndef PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H
#define PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H

#include "AlgorithmsLib_Exports.h"

#include "CandidateClassifier.h"
#include "CandidateProcessertron.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "ParquetReader.h"


using namespace Error;

class CandidateScores;
class TargetDecoyCandidatePair;


namespace NeuralNetDataNamespace {
    const QString SCORES = QStringLiteral("scores");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString IS_TEST = QStringLiteral("isTest");
    const QString PEP_STR = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString TARGET_KEY = QStringLiteral("targetKey");
}

struct NeuralNetData : public ParquetReaderInputBase {

public:

    QVector<double> scores;
    bool isDecoy = false;
    bool isTest = false;
    Charge charge = -1;
    QString targetKey;
    PeptideStringWithMods peptideStringWithMods;

    QMap<QString, QVariant> map() override {

        using namespace NeuralNetDataNamespace;

        return {
                {SCORES, QVariant(qVectorToQByteArray(scores))},
                {IS_DECOY, QVariant(isDecoy)},
                {IS_TEST, QVariant(isTest)},
                {CHARGE, QVariant(charge)},
                {TARGET_KEY, QVariant(targetKey)},
                {PEP_STR, QVariant(peptideStringWithMods)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        ERR_INIT

        using namespace NeuralNetDataNamespace;

        const QStringList keysToCheck = {
                SCORES,
                IS_DECOY,
                IS_TEST,
                CHARGE,
                TARGET_KEY,
                PEP_STR
        };

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        scores = bytesArrayToQVector<double>(dataMap.value(SCORES).toByteArray());
        isDecoy = dataMap.value(IS_DECOY).toBool();
        isTest = dataMap.value(IS_TEST).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR).toString();
        charge = dataMap.value(CHARGE).toInt();
        targetKey = dataMap.value(TARGET_KEY).toString();

        ERR_RETURN
    }
};


class ALGORITHMSLIB_EXPORTS FDRCLassifierNeuralNet {

public:

    FDRCLassifierNeuralNet();
    ~FDRCLassifierNeuralNet();

    Err init(
            int topNMs2Ions,
            int epochs,
            int baggingSize,
            double batchFraction,
            double learningRate,
            const QPair<double, double> &scanTimeMinMax
            );

    Err exec(
            const QVector<ScoredCandidate> &trainingDataTargetsAndDecoys,
            QVector<ScoredCandidate> *scoredCandidatesClassifier
            );

    static QString buildTargetDecoyKey(
            const PeptideStringWithMods &peptideStringWithMods,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            Charge charge
    );

    static QVector<double> buildScoreVector(
            const CandidateScores &candidateScores,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax
    );

    template <typename T>
    static QVector<double> extractScoresFromVecFeatures(
            const QVector<T> &featureVec,
            int theoMzIonsSize
    ) {

        QVector<double> vec(theoMzIonsSize, 0.0);

        for (int i = 0; i < theoMzIonsSize; i++) {

            if (i >= featureVec.size()) {
                break;
            }

            vec[i] = static_cast<double>(featureVec[i]);
        }

        return vec;
    }

    static Err countScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
            );

    static Err countScoreCandidatesByFDR(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePair,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
    );

    static Err outputFDRResults(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
            bool verbose,
            QMap<QString, int> *fdrVsCount
    );

private:

    Err trainClassifier(
            const QMap<QString, ScoredCandidate> &keyVsScoredCandidateCulled,
            QVector<QVector<float>> *allDataVecs,
            QVector<NeuralNetData> *trainingData
            );

    Err trainBaggedNeuralNets(
            const QVector<QVector<QVector<float>>> &trainingDataVecsTranched,
            const QVector<QVector<float>> &ytrainingDataTranched
            );

    Err predictBaggedClassifiers(
            const QVector<QVector<float>> &allDataVecs,
            QVector<float> *meanPredictions
            );

private:

    int m_epochs;
    int m_baggingSize;
    int m_topNMs2Ions;
    double m_batchFraction;
    double m_learningRate;

    const int m_minTopNMs2Ions;
    bool m_isInit;

    QPair<double, double> m_scanTimeMinMax;
    QVector<CandidateClassifier*> m_candidateClassifiers;

};


#endif //PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H
