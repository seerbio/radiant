//
// Created by anichols on 9/2/23.
//

#ifndef PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
#define PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"

using ScoresTargets = QVector<float>;
using ScoresDecoys = QVector<float>;
using ScoresTargetsAndDecoys = QVector<float>;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS ClassifierWeightsManager {

public:

    ClassifierWeightsManager() = default;
    ~ClassifierWeightsManager() = default;

    [[nodiscard]] const QVector<double> &getWeights() const;

    void setWeights(const QVector<double> &weights);

    [[nodiscard]] const QVector<double> &getWeightsBest() const;

    void setWeightsBest(const QVector<double> &weightsBest);

    [[nodiscard]] const QVector<double> &getGuideWeights() const;

    void setGuideWeights(const QVector<double> &guideWeights);

    [[nodiscard]] const QVector<double> &getGuideWeightsBest() const;

    void setGuideWeightsBest(const QVector<double> &guideWeightsBest);

    [[nodiscard]] const QVector<double> &getSelectionWeights() const;

    void setSelectionWeights(const QVector<double> &selectionWeights);

    [[nodiscard]] const QVector<double> &getSelectionWeightsBest() const;

    void setSelectionWeightsBest(const QVector<double> &selectionWeightsBest);

    static Err buildDataClassifier1(
            const QVector<QVector<float>*> &targets,
            const QVector<QVector<float>*> &decoys,
            QVector<QVector<float>> *A,
            QVector<float> *b
            );

    static Err buildDataClassifier2(
            const QVector<QVector<float>*> &targets,
            const QVector<QVector<float>*> &decoys,
            QVector<QVector<float>> *A,
            QVector<float> *b
    );

    static Err fitWeights(
            const QVector<QVector<float>> &matA,
            const QVector<float> &vecB,
            QVector<float> *weights
            );

    static Err applyWeights(
            const QVector<QVector<float>*> &matA,
            const QVector<float> &weights,
            QVector<float> *results
            );

private:

    QVector<double> m_weights;
    QVector<double> m_weightsBest;
    QVector<double> m_guideWeights;
    QVector<double> m_guideWeightsBest;
    QVector<double> m_selectionWeights;
    QVector<double> m_selectionWeightsBest;

};


#endif //PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
