//
// Created by Drucifer on 3/6/2022.
//

#ifndef PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H
#define PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H

#include "AlgorithmsLib_Exports.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "TandemPredictionUtils.h"

using namespace Error;

namespace PeptidePredictionInputNamespace {

    const QString PEPTIDE = QStringLiteral("Peptide");
    const QString CHARGE = QStringLiteral("Charge");
    const QString COLL_ENERGY = QStringLiteral("CollisionEnergy");
    const QString IS_DECOY = QStringLiteral("IsDecoy");

    const QStringList keysToCheck = {
            PEPTIDE,
            CHARGE,
            COLL_ENERGY,
            IS_DECOY
    };
}

struct PeptidePredictionInput : CSVReaderInputBase {

    QString peptideSequence;
    double collisionEnergy = -1.0;
    int charge = -1;
    bool isDecoy = false;

    double normalizedCollisionEnergy = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace PeptidePredictionInputNamespace;

        return {
            {PEPTIDE, peptideSequence},
            {CHARGE, charge},
            {COLL_ENERGY, collisionEnergy},
            {IS_DECOY, isDecoy}
        };
    }

    Err initFromRead(const CSVReaderInputBase &row) override {

        using namespace PeptidePredictionInputNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = CSVReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideSequence = dataMap.value(PEPTIDE).toString();
        collisionEnergy = dataMap.value(COLL_ENERGY).toDouble();
        charge = dataMap.value(CHARGE).toInt();
        isDecoy = dataMap.value(IS_DECOY).toBool();

        ERR_RETURN
    }

};

class ALGORITHMSLIB_EXPORTS TandemFragmentPredictotron {

public:

    using TandemPrediction = QVector<FragmentIon>;

    TandemFragmentPredictotron();
    ~TandemFragmentPredictotron();

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &modelFilePath,
            int charge
            );

    Err  batchPredictTandemSpectra(
            const QVector<PeptidePredictionInput> &predictionInputs,
            QHash<PeptideSequenceChargeKey, TandemPrediction> *tandemPredictions
            );

    void setIntensityThreshold(double intensityThreshold);


private:

    Q_DISABLE_COPY(TandemFragmentPredictotron) class Private;
    const QScopedPointer<Private> d_ptr;

    PythiaParameters m_params;

};


#endif //PYTHIACPP_TANDEMFRAGMENTPREDICTOTRON_H
