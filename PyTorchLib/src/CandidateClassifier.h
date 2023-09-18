//
// Created by anichols on 9/17/23.
//

#ifndef PYTHIADIACPP_CANDIDATECLASSIFIER_H
#define PYTHIADIACPP_CANDIDATECLASSIFIER_H

#include "PyTorchLib_Exports.h"


#include <QVector>


class PYTORCHLIB_EXPORTS CandidateClassifier {

public:

    CandidateClassifier() = default;
    ~CandidateClassifier() = default;

    bool trainCandidateClassifier(
            const QVector<QVector<double>> &xData,
            const QVector<double> &yData
            );


};


#endif //PYTHIADIACPP_CANDIDATECLASSIFIER_H
