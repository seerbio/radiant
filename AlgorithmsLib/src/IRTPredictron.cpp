//
// Created by anichols on 7/22/23.
//

#include "IRTPredictron.h"

#include "AminoAcids.h"
#include "CalibrationReader.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "NeuralNetModel.h"
#include "TandemPredictionUtils.h"

#include <Eigen/Dense>


class Q_DECL_HIDDEN IRTPredictron::Private
{
public:

    Private();
    ~Private();

    Err init(const QString &modelFilePath);

    Err batchPredictIRT(
            const QStringList &peptideStringWithModsList,
            QVector<float> *rawPredictionResults
            );

private:

    NeuralNetModel *m_model;
    bool m_isInit;
    const int m_maxSequenceLength;

};

IRTPredictron::Private::Private()
: m_isInit(false)
, m_model(nullptr)
, m_maxSequenceLength(40)
{}

IRTPredictron::Private::~Private() {
    delete m_model;
}

Err IRTPredictron::Private::init(const QString &modelFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(modelFilePath); ree

    m_model = new NeuralNetModel();
    e = m_model->init(modelFilePath); ree;

    m_isInit = true;

    ERR_RETURN
}

namespace{

    bool checkIfPeptideLengthsAreValid(
            const QStringList &peptideStringWithModsList,
            int maxPeptideLength
            ) {

        const auto checkingLogic = [&](const QString &s){
            return s.size() > maxPeptideLength;
        };

        return std::none_of(peptideStringWithModsList.begin(),  peptideStringWithModsList.end(), checkingLogic);
    }

    Err convertPeptideSequencesToOneHotMatrices(
            const QStringList &peptideStringWithModsList,
            int maxPeptideSize,
            QVector<Eigen::MatrixXd> *mats
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithModsList); ree

        const int aminoAcidCountPlusNull = 20 + 1;

        for (const QString &peptideStringWithMods : peptideStringWithModsList) {

            const QVector<int> seqToNumericalVec = TandemPredictionUtils::peptideSequenceToArray(
                    peptideStringWithMods,
                    maxPeptideSize
                    );

            Eigen::MatrixXd mat = EigenUtils::vectorToOneHotMatrix(
                    seqToNumericalVec,
                    aminoAcidCountPlusNull
                    );

            mats->push_back(mat);
        }

        ERR_RETURN
    }

}//namespace
Err IRTPredictron::Private::batchPredictIRT(
        const QStringList &peptideStringWithModsList,
        QVector<float> *rawPredictionResults
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(peptideStringWithModsList); ree

    const bool allSequenceLengthsAreValid = checkIfPeptideLengthsAreValid(
            peptideStringWithModsList,
            m_maxSequenceLength
            );
    e = ErrorUtils::isTrue(allSequenceLengthsAreValid); ree

    QVector<Eigen::MatrixXd> mats;
    e = convertPeptideSequencesToOneHotMatrices(
            peptideStringWithModsList,
            m_maxSequenceLength,
            &mats
            ); ree

    const QVector<QVector<float>> results = m_model->batchPredict(mats);

    std::transform(
            results.begin(),
            results.end(),
            std::back_inserter(*rawPredictionResults),
            [](const QVector<float> &r){return r.front();}
            );

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


IRTPredictron::IRTPredictron() : d_ptr(new Private()) {}

IRTPredictron::~IRTPredictron() {}

Err IRTPredictron::init(const QString &modelFilePath) {
    ERR_INIT
    e = d_ptr->init(modelFilePath); ree
    ERR_RETURN
}

Err IRTPredictron::batchPredictIRT(
        const QStringList &peptideStringWithModsList,
        QVector<float> *rawPredictionResults
        ) {

    ERR_INIT
    e = d_ptr->batchPredictIRT(peptideStringWithModsList, rawPredictionResults); ree
    ERR_RETURN
}

//Err IRTPredictron::buildNearestNeighborsIRTData(
//        const QString &iRTRecalibrationFilePath,
//        QVector<QPair<double, Coors>> *nnInputData
//        ) {
//
//    ERR_INIT
//
//    e = ErrorUtils::fileExists(iRTRecalibrationFilePath); ree
//
//    QVector<IRTReCalibrationRow> iRTReCalibrationReaderRows;
//    e  = CSVReader::read(
//            iRTRecalibrationFilePath,
//            &iRTReCalibrationReaderRows
//    ); ree;
//
//    for (const IRTReCalibrationRow &row : iRTReCalibrationReaderRows) {
//        nnInputData->push_back({row.scanTime, {row.iRT, 0.0}});
//    }
//
//    ERR_RETURN
//
//}
