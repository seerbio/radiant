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
#include "PeptideStringWithMods.h"


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

/**
 * See ParquetReaderInputBase for documentation
 */
struct NeuralNetData : public ParquetReaderInputBase {

public:

    QVector<double> scores;
    double discScore;
    bool isDecoy = false;
    int charge = -1;
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
        peptideStringWithMods = PeptideStringWithMods(dataMap.value(PEP_STR).toString());
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

    /**
    * @brief Initializes the FDR classifier with neural networks.
    *
    * This function initializes the FDR classifier with neural networks by setting various parameters.
    *
    * @param epochs An integer representing the number of training epochs.
    * @param baggingSize An integer representing the size of the bagging ensemble.
    * @param batchSize An integer representing the batch size for training.
    * @param learningRate A double representing the learning rate for training.
    * @return An Err enum indicating the success or failure of the operation.
    */
    Err init(
            int epochs,
            int baggingSize,
            int batchSize,
            double learningRate,
            int threadCount
            );

    /**
    * @brief Executes the FDR classifier using neural networks.
    *
    * This function performs the execution of the FDR classifier using neural networks. It includes
    * training the classifier with input data (xData, yData) and generating mean predictions.
    *
    * @param xData A QVector of QVector<float> representing the input data for training and prediction.
    * @param yData A QVector<float> representing the target values for training.
    * @param seed An integer representing the seed for random number generation.
    * @param meanPredictions Output parameter, a QVector<float>* representing the mean predictions.
    * @return An Err enum indicating the success or failure of the operation.
    */
    Err exec(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int seed,
            int verbosity,
            QVector<float> *meanPredictions
    );

    Err trainClassifier(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int seed,
            int verbosity
            );

    Err predictBaggedClassifiers(
            const QVector<QVector<float>> &allDataVecs,
            QVector<float> *meanPredictions
            );


    /**
    * @brief Counts the number of target candidates below a specified q-value threshold.
    *
    * This function takes a QVector of CandidateScores and counts the number of target candidates
    * with a q-value below the specified threshold.
    *
    * @param candidateScores A QVector of CandidateScores representing the candidate scores.
    * @param qValueThreshold A double representing the q-value threshold.
    * @param targetCountBelowFDRThreshold Output parameter, an integer pointer storing the count of target candidates below the threshold.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err countScoreCandidatesByFDR(
            const QVector<CandidateScores> &targetDecoyCandidatePair,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
    );

    static Err countScoreCandidatesByFDR(
            QVector<CandidateScores*> &targetDecoyCandidatePair,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
    );

    /**
    * @brief Outputs FDR results based on specified FDR thresholds.
    *
    * This function outputs FDR results for different FDR thresholds using the given candidate scores.
    *
    * @param candidateScores A QVector of CandidateScores representing the candidate scores.
    * @param verbose A boolean indicating whether to display verbose output.
    * @param fdrVsCount Output parameter, a QMap<QString, int>* mapping FDR percentages to counts of candidates below each threshold.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err outputFDRResults(
            const QVector<CandidateScores> &candidateScores,
            bool verbose,
            QMap<int, int> *fdrVsCount
    );

    static Err outputFDRResults(
            QVector<CandidateScores*> &candidateScores,
            int verbose,
            QMap<int, int> *fdrVsCount
    );

    static Err outPutFDRCounts(
        const QMap<int, int> &fdrVsCount,
        QString *outputString
        );


private:

    Err trainBaggedNeuralNets(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int seed,
            int verbosity
            );



private:

    int m_epochs;
    int m_baggingSize;
    int m_batchSize;
    double m_learningRate;
    int m_threadCount;

    bool m_isInit;

    QVector<CandidateClassifier*> m_candidateClassifiers;

};


#endif //PYTHIADIACPP_FDRCLASSIFIERNEURALNET_H
