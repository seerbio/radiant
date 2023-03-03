//
// Created by Drucifer on 3/6/2022.
//

#ifndef PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H
#define PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "TandemPredictionUtils.h"

using namespace Error;

struct PeptidePredictionInput {
    QString peptideSequence;
    double collisionEnergy = -1.0;
    int charge = -1;

    double normalizedCollisionEnergy = -1.0;
};

class ALGORITHMSLIB_EXPORTS TandemFragmentPredictotron {

public:

    using TandemPrediction = QVector<FragmentIon>;

    TandemFragmentPredictotron();
    ~TandemFragmentPredictotron();

    Err init(
            const QString &modelFilePath,
            int charge
            );

    Err  batchPredictTandemSpectra(
            const QVector<PeptidePredictionInput> &predictionInputs,
            QHash<PeptideSequenceChargeCollisionEnergyKey, TandemPrediction> *tandemPredictions
            );

    void setIntensityThreshold(double intensityThreshold);


private:

    Q_DISABLE_COPY(TandemFragmentPredictotron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H
