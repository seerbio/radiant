//
// Created by Drucifer on 3/6/2022.
//

#include "TandemFragmentPredictotron.h"

#include "AminoAcids.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "NeuralNetModel.h"

#include <Eigen/Dense>


class Q_DECL_HIDDEN TandemFragmentPredictotron::Private
{
public:

    Private();
    ~Private();

    Err init(
            const QString &modelFilePath,
            int charge
            );

    Err  batchPredictTandemSpectra(
            const QVector<PeptidePredictionInput> &tandemPredictionInputs,
            QHash<PeptideSequenceChargeKey, TandemPrediction> *result
            );

    void setIntensityThreshold(double intensityThreshold);


private:

    NeuralNetModel* m_model;
    int m_charge;
    bool m_isInit;
    QStringList m_ionLabels;
    double m_intensityThreshold;

};

TandemFragmentPredictotron::Private::Private()
: m_isInit(false)
, m_model(nullptr)
, m_charge(-1)
, m_ionLabels({})
, m_intensityThreshold(0.01) {}

TandemFragmentPredictotron::Private::~Private() {
    if(m_isInit){
        delete m_model;
    }
}

Err TandemFragmentPredictotron::Private::init(
        const QString &modelFilePath,
        int charge
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(modelFilePath); ree;
    e = ErrorUtils::fileExists(modelFilePath); ree;

    auto *model = new NeuralNetModel();
    e = model->init(modelFilePath); ree;
    m_model = model;

    const int minCharge = 1;
    const int maxCharge = 4;

    m_charge = charge;
    e = ErrorUtils::isWithinRange(
            m_charge,
            minCharge,
            maxCharge,
            ErrorUtilsParam::IncludeThreshold
            ); ree;

    m_ionLabels = TandemPredictionUtils::buildIonLabels(m_charge);
    e = ErrorUtils::isNotEmpty(m_ionLabels); ree;

    m_isInit = true;

    ERR_RETURN
}

namespace {

    // These methods are in this namespace so Eigen is not exposed to header if they were in TandemPredictionUtils.

    bool checkIfPeptideLengthsAreValid(
            const QVector<PeptidePredictionInput> &predictionInputs,
            int maxPeptideLength
            ) {

        const auto checkingLogic = [&](const PeptidePredictionInput &s){
            return s.peptideSequence.size() > maxPeptideLength;
        };

        return std::none_of(predictionInputs.begin(),  predictionInputs.end(), checkingLogic);
    }

    Err convertPeptideSequencesToOneHotMatrices(
            const QVector<PeptidePredictionInput> &predictionInputs,
            int maxPeptideSize,
            QVector<Eigen::MatrixXd> *mats
            ) {

        ERR_INIT

        const int aminoAcidCountPlusNull = 20 + 1;
        const double normalizerNCE = 100.0;

        for (const PeptidePredictionInput &predictionInput : predictionInputs) {

            const QVector<int> seqToNumericalVec
                        = TandemPredictionUtils::peptideSequenceToArray(predictionInput.peptideSequence, maxPeptideSize);

            Eigen::MatrixXd mat = EigenUtils::vectorToOneHotMatrix(seqToNumericalVec, aminoAcidCountPlusNull);

            mat = mat.array() * (predictionInput.normalizedCollisionEnergy / normalizerNCE);
            mats->push_back(mat);
        }

        ERR_RETURN
    }

    QVector<QVector<float>> normalizePredictions(const QVector<QVector<float>> &predictions) {

        QVector<QVector<float>> returnVecs;
        for (QVector<float> vec : predictions) {
            //TODO replace w/ EigenUtils::convertQvectorToEigenVector
            Eigen::VectorX<float> eVec = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(vec.data(), vec.size());
            EigenUtils::normalizeVector<Eigen::VectorX<float>, float>(&eVec);
            EigenUtils::thresholdVector(static_cast<float>(1e-2), &eVec);

            std::vector<float> vecReturn(eVec.data(), eVec.data() + eVec.size());
            returnVecs.push_back(QVector<float>::fromStdVector(vecReturn));
        }

        return returnVecs;
    }

    Err zipFragmentIonList(
            const QVector<float> &intensityVals,
            const QStringList &ionLabels,
            QVector<FragmentIon> *frags
            ) {

            ERR_INIT

            e = ErrorUtils::isEqual(intensityVals.size(), ionLabels.size()); ree;

            for (int i = 0; i < intensityVals.size(); i++) {

                const float intensity = MathUtils::pRound(intensityVals.at(i), 2);
                const QString &ionLabel = ionLabels.at(i);

                FragmentIon fragmentIon(intensity, ionLabel);
                frags->push_back(fragmentIon);
            }

            ERR_RETURN
    }

    QString buildPeptideLookupKey(const PeptidePredictionInput &ppi) {
        return ppi.peptideSequence + S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP + QString::number(ppi.charge);
    }


}//NAMESPACE
Err  TandemFragmentPredictotron::Private::batchPredictTandemSpectra(
        const QVector<PeptidePredictionInput> &predictionInputs,
        QHash<PeptideSequenceChargeKey, TandemPrediction> *tandemPredictions
        ) {

    ERR_INIT

    const int maxSequenceLength = TandemPredictionUtils::getMaxLengthForPredictionInputByCharge(m_charge);

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(predictionInputs); ree;
    e = ErrorUtils::isNotEmpty(m_ionLabels); ree;

    const bool allSequenceLengthsAreValid = checkIfPeptideLengthsAreValid(predictionInputs, maxSequenceLength);
    e = ErrorUtils::isTrue(allSequenceLengthsAreValid); ree;

    QVector<Eigen::MatrixXd> mats;
    e = convertPeptideSequencesToOneHotMatrices(predictionInputs, maxSequenceLength, &mats); ree;

    const QVector<QVector<float>> rawPredictionResults = m_model->batchPredict(mats);
    const QVector<QVector<float>> rawPredictionResultsCleaned = normalizePredictions(rawPredictionResults);

    for (int i = 0; i < rawPredictionResultsCleaned.size(); i++) {

        const QVector<float> &intensityVals = rawPredictionResultsCleaned.at(i);

        const PeptidePredictionInput &ppi = predictionInputs.at(i);

        const PeptideSequenceChargeKey peptideLookupKey = buildPeptideLookupKey(ppi);

        TandemPrediction frags;
        e = zipFragmentIonList(
                intensityVals,
                m_ionLabels,
                &frags
                ); ree;

        TandemPredictionUtils::filterFragmentIonVectorByIntensity(m_intensityThreshold, &frags);
        TandemPredictionUtils::sortPredictionByIonLabel(&frags);
        tandemPredictions->insert(peptideLookupKey, frags); ree;
    }

    ERR_RETURN
}

void TandemFragmentPredictotron::Private::setIntensityThreshold(double intensityThreshold) {
    m_intensityThreshold = intensityThreshold;
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


TandemFragmentPredictotron::TandemFragmentPredictotron() : d_ptr(new Private()) {}

TandemFragmentPredictotron::~TandemFragmentPredictotron() {}

Err TandemFragmentPredictotron::init(
        const QString &modelFilePath,
        int charge
        ) {
    ERR_INIT
    e = d_ptr->init(modelFilePath, charge); ree;
    ERR_RETURN
}

Err TandemFragmentPredictotron::batchPredictTandemSpectra(
        const QVector<PeptidePredictionInput> &tandemPredictionInputs,
        QHash<PeptideSequenceChargeKey, TandemPrediction> *result
        ) {
   ERR_INIT
   e = d_ptr->batchPredictTandemSpectra(tandemPredictionInputs, result); ree;
   ERR_RETURN
}

void TandemFragmentPredictotron::setIntensityThreshold(double intensityThreshold) {
    d_ptr->setIntensityThreshold(intensityThreshold);
}
