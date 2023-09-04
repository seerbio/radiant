//
// Created by anichols on 9/2/23.
//

#ifndef PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
#define PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"

using namespace Error;

class ALGORITHMSLIB_EXPORTS ClassifierWeightsManager {

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
            const QVector<QVector<double>> &targets,
            const QVector<QVector<double>> &decoys,
            QVector<QVector<double>> *A,
            QVector<double> *b
            );

    static Err fitWeights(
            const QVector<QVector<double>> &matA,
            const QVector<double> &vecB,
            QVector<double> *weights
            );

    static Err applyWeights(
            const QVector<QVector<double>> &matA,
            const QVector<double> &weights,
            QVector<double> *results
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
