//
// Created by anichols on 9/17/23.
//

#ifndef PYTHIADIACPP_CANDIDATECLASSIFIER_H
#define PYTHIADIACPP_CANDIDATECLASSIFIER_H

#include "PyTorchLib_Exports.h"

#include <QScopedPointer>
#include <QVector>


class PYTORCHLIB_EXPORTS CandidateClassifier {

public:

    CandidateClassifier();
    ~CandidateClassifier();

    bool trainCandidateClassifier(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int epochsMax,
            double batchFraction,
            double learningRate
            );

private:

    Q_DISABLE_COPY(CandidateClassifier) class Private;
    const QScopedPointer<Private> d_ptr;


};


#endif //PYTHIADIACPP_CANDIDATECLASSIFIER_H
