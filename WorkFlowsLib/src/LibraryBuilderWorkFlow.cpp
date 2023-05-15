//
// Created by anichols on 2/17/23.
//

#include "LibraryBuilderWorkFlow.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "MolecularFormula.h"
#include "PythiaParameterReader.h"
#include "TandemFragmentPredictotron.h"
#include "FragLibReader.h"
#include "TandemPredictionUtils.h"

#include <QtConcurrent/QtConcurrent>


LibraryBuilderWorkFlow::~LibraryBuilderWorkFlow() {

    for (int i = 1; i < 5; i++) {
        delete m_tandemPredictionModels[i];
    }

}

Err LibraryBuilderWorkFlow::init(
        const PythiaParameters &pythiaParameters,
        const QString &modelCharge1,
        const QString &modelCharge2,
        const QString &modelCharge3,
        const QString &modelCharge4
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;

    m_modelFilePaths.insert(1, modelCharge1);
    m_modelFilePaths.insert(2, modelCharge2);
    m_modelFilePaths.insert(3, modelCharge3);
    m_modelFilePaths.insert(4, modelCharge4);

    for (auto it = m_modelFilePaths.begin(); it != m_modelFilePaths.end(); it++) {

        const int charge = it.key();
        const QString &modelFilePath = it.value();

        auto *model = new TandemFragmentPredictotron();
        e = model->init(
                m_pythiaParameters,
                m_modelFilePaths.value(charge),
                charge
                ); ree;

        m_tandemPredictionModels.insert(charge, model);
    }

    ERR_RETURN
}


namespace {

    Err buildTandemPredictionInputs(
            const AminoAcids &aminoAcids,
            QVector<PeptidePredictionInput> *peptidePredictionInputs
            ) {

        ERR_INIT

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

    Err writePredictionsToCSV(
            const QHash<PeptideSequenceChargeKey,
                    TandemFragmentPredictotron::TandemPrediction> &tandemPredictionsAllCharges,
            const QString &outputFilePath
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(outputFilePath); ree;
        e = ErrorUtils::isNotEmpty(tandemPredictionsAllCharges); ree;

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

                const PeptideSequenceChargeKey &key = itt.key();
                TandemFragmentPredictotron::TandemPrediction pred = itt.value();

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

    QVector<FragLibReaderRow> buildFragLibReaderRows(
            const QHash<PeptideSequenceChargeKey,TandemFragmentPredictotron::TandemPrediction> &tandemPredictionsAllCharges,
            const QHash<PeptideSequenceChargeKey, bool> &peptideSeqChargeKeyVsIsDecoy,
            const AminoAcids &aminoAcids
            ) {

        QVector<FragLibReaderRow> tlrs;
        for (auto it = tandemPredictionsAllCharges.begin(); it != tandemPredictionsAllCharges.end(); it++) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            const TandemFragmentPredictotron::TandemPrediction &tandemPrediction = it.value();

            QStringList ionLabels;
            QVector<double> mzVals;
            QVector<double> intensities;
            for (const FragmentIon &fi : tandemPrediction) {
                mzVals.push_back(fi.mz);
                ionLabels.push_back(fi.ionLabel);
                intensities.push_back(fi.intensity);
            }

            const PeptideString peptideString = peptideSequenceChargeKey.split(
                    S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                    QString::SkipEmptyParts
                    ).front();

            FragLibReaderRow row;
            row.peptideSequenceChargeKey = peptideSequenceChargeKey;
            row.mzVals = mzVals;
            row.mass = BiophysicalCalcs::calculatePeptideMass(peptideString, aminoAcids);
            row.intensityVals = intensities;
            row.isDecoy = peptideSeqChargeKeyVsIsDecoy.value(row.peptideSequenceChargeKey);

            tlrs.push_back(row);
        }

        const auto sortLogicMassAsc = [](
                const FragLibReaderRow &l,
                const FragLibReaderRow &r
                ){
            return l.mass < r.mass;
        };

        std::sort(tlrs.begin(), tlrs.end(), sortLogicMassAsc);

        return tlrs;
    }

    Err writePredictionsToParquet(
            const QHash<PeptideSequenceChargeKey, TandemFragmentPredictotron::TandemPrediction> &tandemPredictionsAllCharges,
            const QString &outputFilePath,
            const QHash<PeptideSequenceChargeKey, bool> &peptideSequenceChargeKeyVsIsBool,
            const AminoAcids &aminoAcids
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(outputFilePath); ree;
        e = ErrorUtils::isNotEmpty(tandemPredictionsAllCharges); ree;

        const QVector<FragLibReaderRow> writeRows = buildFragLibReaderRows(
                tandemPredictionsAllCharges,
                peptideSequenceChargeKeyVsIsBool,
                aminoAcids
                );

        e = ParquetReader::write(
                writeRows,
                outputFilePath
                ); ree;

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
    e = CSVReader::read(
            peptidesCSVFilePath,
            &peptidePredictionInputs
            ); ree;

    e = buildPeptideSequenceChargeKeyVsIsDecoy(peptidePredictionInputs); ree;

    e = buildTandemPredictionInputs(
            m_pythiaParameters.aminoAcids,
            &peptidePredictionInputs
            ); ree;

    const QMap<Charge, QVector<PeptidePredictionInput>> sortedPeptidePredictionInputs
            = sortPeptidePredictionInputs(peptidePredictionInputs);

    QHash<PeptideSequenceChargeKey,
            TandemFragmentPredictotron::TandemPrediction> tandemPredictionsAllCharges;

    for (auto it = sortedPeptidePredictionInputs.begin(); it != sortedPeptidePredictionInputs.end(); it++) {

        const Charge charge = it.key();
        const QVector<PeptidePredictionInput> &ppis = it.value();

        qDebug() << "Prediction charge:" << charge;

        TandemFragmentPredictotron *model = m_tandemPredictionModels.value(charge);

        QHash<PeptideSequenceChargeKey,
              TandemFragmentPredictotron::TandemPrediction> tandemPredictions;

        e = model->batchPredictTandemSpectra(ppis, &tandemPredictions); ree;

        for (auto itt = tandemPredictions.begin(); itt != tandemPredictions.end(); itt++) {

            const PeptideSequenceChargeKey &key = itt.key();
            const TandemFragmentPredictotron::TandemPrediction &pred = itt.value();

            tandemPredictionsAllCharges.insert(key, pred);
        }

    }

    *returnFilePath = peptidesCSVFilePath + S_GLOBAL_SETTINGS.DOT_FRAGLIB;

    e = writePredictionsToParquet(
            tandemPredictionsAllCharges,
            *returnFilePath,
            m_peptideSequenceChargeKeyVsIsDecoy,
            m_pythiaParameters.aminoAcids
            ); ree;

    ERR_RETURN
}

Err LibraryBuilderWorkFlow::buildPeptideSequenceChargeKeyVsIsDecoy(
        const QVector<PeptidePredictionInput> &peptidePredictionInputs) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(peptidePredictionInputs); ree;

    for (const PeptidePredictionInput &ppi : peptidePredictionInputs) {

        const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
                ppi.peptideSequence,
                ppi.charge
        );

        m_peptideSequenceChargeKeyVsIsDecoy.insert(peptideSequenceChargeKey, ppi.isDecoy);

    }

    ERR_RETURN
}
