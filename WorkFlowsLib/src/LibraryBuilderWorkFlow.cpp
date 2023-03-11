//
// Created by anichols on 2/17/23.
//

#include "LibraryBuilderWorkFlow.h"

#include "BiophysicalCalcs.h"
#include "FragLibraryTron.h"
#include "MolecularFormula.h"
#include "TandemFragmentPredictotron.h"
#include "TandemPredictionUtils.h"

#include <QtConcurrent/QtConcurrent>


LibraryBuilderWorkFlow::~LibraryBuilderWorkFlow() {

    for (int i = 1; i < 5; i++) {
        delete m_tandemPredictionModels[i];
    }

}

Err LibraryBuilderWorkFlow::init(
        const QString &modelCharge1,
        const QString &modelCharge2,
        const QString &modelCharge3,
        const QString &modelCharge4
        ) {

    ERR_INIT

    m_modelFilePaths.insert(1, modelCharge1);
    m_modelFilePaths.insert(2, modelCharge2);
    m_modelFilePaths.insert(3, modelCharge3);
    m_modelFilePaths.insert(4, modelCharge4);

    for (auto it = m_modelFilePaths.begin(); it != m_modelFilePaths.end(); it++) {

        const int charge = it.key();
        const QString &modelFilePath = it.value();

        auto *model = new TandemFragmentPredictotron();
        e = model->init(
                m_modelFilePaths.value(charge),
                charge
                ); ree;

        m_tandemPredictionModels.insert(charge, model);
    }

    ERR_RETURN
}


namespace {

    Err readCSV(
            const QString &peptidesCSVFilePath,
            QVector<PeptidePredictionInput> *peptidePredictionInputs
            ) {

        ERR_INIT

        QFile file(peptidesCSVFilePath);
        if (!file.open(QIODevice::ReadOnly)) {
            qDebug() << file.errorString();
            rrr(eFileError);
        }

        const QString expectedHeader = "Peptide,Charge,CollisionEnergy\n";

        QString header;

        while (!file.atEnd()) {

            if (header.isEmpty()) {
                header = file.readLine();

                e = ErrorUtils::isEqualString(header, expectedHeader); ree;
                continue;
            }

            const QString &row = file.readLine().replace(S_GLOBAL_SETTINGS.NEWLINE, "");
            const QStringList rowSplit = row.split(S_GLOBAL_SETTINGS.COMMA);

            PeptidePredictionInput ppi;

            ppi.peptideSequence = rowSplit.at(0);
            e = ErrorUtils::toInt(rowSplit.at(1), &ppi.charge); ree;
            e = ErrorUtils::toDouble(rowSplit.at(2), &ppi.collisionEnergy); ree;

            if (!TandemPredictionUtils::isValidPeptideForPrediction(ppi.peptideSequence, ppi.charge)) {
                continue;
            }

            peptidePredictionInputs->push_back(ppi);
        }

        ERR_RETURN
    }

    Err buildTandemPredictionInputs(QVector<PeptidePredictionInput> *peptidePredictionInputs) {

        ERR_INIT

        AminoAcids aminoAcids;
        aminoAcids.addFixedModification('C', MolecularFormulas::carbamidomethylFormula);

        for (int i = 0; i < peptidePredictionInputs->size(); i++) {

            PeptidePredictionInput &ppi = (*peptidePredictionInputs)[i];

            e = ErrorUtils::isTrue(ppi.charge > 0); ree;

            const double mass = BiophysicalCalcs::calculatePeptideMass(ppi.peptideSequence, aminoAcids);
            const double mz = BiophysicalCalcs::calculateThomsonFromMass(mass, ppi.charge);

            ppi.normalizedCollisionEnergy = TandemPredictionUtils::calculateNormalizedCollisionEnergy(
                    mz,
                    ppi.charge,
                    static_cast<int>(ppi.collisionEnergy)
                    );
        }

        ERR_RETURN
    }

    QMap<Charge, QVector<PeptidePredictionInput>> sortPeptidePredictionInputs(
            const QVector<PeptidePredictionInput> &peptidePredictionInputs
            ) {

        QMap<Charge, QVector<PeptidePredictionInput>> sortedPeptidePredictionInputs;

        for (const PeptidePredictionInput &ppi : peptidePredictionInputs) {

            sortedPeptidePredictionInputs[ppi.charge].push_back(ppi);
        }

        return sortedPeptidePredictionInputs;
    };

    void takeTop8IntensityFrags(TandemFragmentPredictotron::TandemPrediction *pred) {

        const auto logic = [](const FragmentIon &l, const FragmentIon &r){
            return l.intensity < r.intensity;
        };

        std::sort(pred->rbegin(), pred->rend(), logic);

        pred->resize(std::min(8, pred->size()));
    }

    Err writePredictionsToCSV(
            const QHash<PeptideSequenceChargeCollisionEnergyKey,
            TandemFragmentPredictotron::TandemPrediction> &tandemPredictionsAllCharges,
            const QString &outputFilePath
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(outputFilePath); ree;

        const QStringList rowHeader = {
                QStringLiteral("Peptide"),
                QStringLiteral("IonLabel"),
                QStringLiteral("Intensity")
        };

        QFile file(outputFilePath);
        if (file.open(QIODevice::ReadWrite)) {

            QTextStream stream(&file);
            stream << rowHeader.join(S_GLOBAL_SETTINGS.COMMA) + S_GLOBAL_SETTINGS.NEWLINE;

            for (auto itt = tandemPredictionsAllCharges.begin(); itt != tandemPredictionsAllCharges.end(); itt++) {

                const PeptideSequenceChargeCollisionEnergyKey &key = itt.key();
                TandemFragmentPredictotron::TandemPrediction pred = itt.value();

                takeTop8IntensityFrags(&pred);

                for (const FragmentIon &fli : pred) {

                    const QStringList row = {
                            key,
                            fli.ionLabel,
                            QString::number(fli.intensity)};

                    stream << row.join(S_GLOBAL_SETTINGS.COMMA) + S_GLOBAL_SETTINGS.NEWLINE;
                }
            }
        }

        qDebug() << "FragLib file written to" << outputFilePath;
        file.close();
        ERR_RETURN
    }

}//namespace
Err LibraryBuilderWorkFlow::exec(
        const QString &peptidesCSVFilePath,
        QString *returnFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_tandemPredictionModels); ree;

    QVector<PeptidePredictionInput> peptidePredictionInputs;
    e = readCSV(
            peptidesCSVFilePath,
            &peptidePredictionInputs
            ); ree;

    e = buildTandemPredictionInputs(&peptidePredictionInputs); ree;

    const QMap<Charge, QVector<PeptidePredictionInput>> sortedPeptidePredictionInputs
            = sortPeptidePredictionInputs(peptidePredictionInputs);

    QHash<PeptideSequenceChargeCollisionEnergyKey,
            TandemFragmentPredictotron::TandemPrediction> tandemPredictionsAllCharges;

    for (auto it = sortedPeptidePredictionInputs.begin(); it != sortedPeptidePredictionInputs.end(); it++) {

        const Charge charge = it.key();
        const QVector<PeptidePredictionInput> &ppis = it.value();

        qDebug() << "Prediction charge:" << charge;

        TandemFragmentPredictotron *model = m_tandemPredictionModels.value(charge);

        QHash<PeptideSequenceChargeCollisionEnergyKey,
              TandemFragmentPredictotron::TandemPrediction> tandemPredictions;

        e = model->batchPredictTandemSpectra(ppis, &tandemPredictions); ree;

        for (auto itt = tandemPredictions.begin(); itt != tandemPredictions.end(); itt++) {

            const PeptideSequenceChargeCollisionEnergyKey &key = itt.key();
            const TandemFragmentPredictotron::TandemPrediction &pred = itt.value();

            tandemPredictionsAllCharges.insert(key, pred);
        }

    }

    *returnFilePath = peptidesCSVFilePath + S_GLOBAL_SETTINGS.DOT_FRAGLIB;

    e = writePredictionsToCSV(
            tandemPredictionsAllCharges,
            *returnFilePath
            ); ree;

    ERR_RETURN
}




