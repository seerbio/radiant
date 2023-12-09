//
// Created by anichols on 12/8/23.
//

#ifndef PYTHIADIACPP_FEATUREFINDERHILLRTREE_H
#define PYTHIADIACPP_FEATUREFINDERHILLRTREE_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class FeatureFinderHill;

class ALGORITHMSFFLIB_EXPORTS FeatureFinderHillRTree {

public:

    FeatureFinderHillRTree();
    ~FeatureFinderHillRTree();

    Err init(const QVector<FeatureFinderHill*> &featureFinderHills);

    Err getHills(
            double mzMin,
            double mzMax,
            FrameIndex frameIndexMin,
            FrameIndex frameIndexMax,
            QVector<FeatureFinderHill*> *featureFinderHills
    );

private:

    Q_DISABLE_COPY(FeatureFinderHillRTree) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_FEATUREFINDERHILLRTREE_H
