//
// Created by anichols on 9/26/23.
//

#ifndef PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H
#define PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateClassifier.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "ParquetReader.h"


using namespace Error;

class CandidateScores;
class TargetDecoyCandidatePair;


namespace NeuralNetDataNamespace {
    const QString SCORES = QStringLiteral("scores");
    const QString DISC_SCORE = QStringLiteral("discScore");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString PEP_STR = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString TARGET_KEY = QStringLiteral("targetKey");
}

struct NeuralNetData : public ParquetReaderInputBase {

public:

    QVector<double> scores;
    double discScore;
    bool isDecoy = false;
    Charge charge = -1;
    QString targetKey;
    PeptideStringWithMods peptideStringWithMods;

    QMap<QString, QVariant> map() override {

        using namespace NeuralNetDataNamespace;

        return {
                {SCORES, QVariant(qVectorToQByteArray(scores))},
                {DISC_SCORE, QVariant(discScore)},
                {IS_DECOY, QVariant(isDecoy)},
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
                DISC_SCORE,
                IS_DECOY,
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
        peptideStringWithMods = dataMap.value(PEP_STR).toString();
        charge = dataMap.value(CHARGE).toInt();
        targetKey = dataMap.value(TARGET_KEY).toString();
        discScore = dataMap.value(DISC_SCORE).toDouble();

        ERR_RETURN
    }
};


class ALGORITHMSFFLIB_EXPORTS FDRCLassifierNeuralNet {

public:

    FDRCLassifierNeuralNet();
    ~FDRCLassifierNeuralNet();

    Err init(
            int epochs,
            int baggingSize,
            int batchSize,
            double learningRate
            );

    Err exec(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            QVector<float> *meanPredictions
    );


    static QString buildTargetDecoyKey(
            const PeptideStringWithMods &peptideStringWithMods,
            const MzTargetKey &uniqueMsInfoScanKey,
            Charge charge
    );

    static Err buildScoreVector(
            const CandidateScores &candidateScores,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax,
            QVector<double> *scoreVec
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
            const QVector<CandidateScores> &targetDecoyCandidatePair,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
    );

    static Err outputFDRResults(
            const QVector<CandidateScores> &candidateScores,
            bool verbose,
            QMap<QString, int> *fdrVsCount
    );

    static Err filterScoreCandidatesByFDR(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
            double qValueThreshold,
            QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairsFDRThresholded
            );

    static Err filterScoreCandidatesByFDR(
            const QVector<CandidateScores> &targetDecoyCandidatePairs,
            double qValueThreshold,
            QVector<CandidateScores> *targetDecoyCandidatePairsFDRThresholded
    );

private:

    Err trainClassifier(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData
            );

    Err trainBaggedNeuralNets(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData
            );

    Err predictBaggedClassifiers(
            const QVector<QVector<float>> &allDataVecs,
            QVector<float> *meanPredictions
            );

private:

    int m_epochs;
    int m_baggingSize;
    int m_topNMs2Ions;
    int m_batchSize;
    double m_learningRate;

    const int m_minTopNMs2Ions;
    bool m_isInit;

    QPair<double, double> m_scanTimeMinMax;
    QVector<CandidateClassifier*> m_candidateClassifiers;

};


#endif //PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H
