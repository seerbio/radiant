//
// Created by anichols on 4/10/23.
//

#include "FragLibraryTronDIA.h"

#include "ErrorUtils.h"
#include "TandemLibraryReader.h"
#include "TandemPredictionUtils.h"

#include <QtConcurrent/QtConcurrent>


Err FragLibraryTronDIA::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    m_params = pythiaParameters;

    buildChargeVsIonLabels();

    e = readFragLibFile(fragLibFilePath); ree;

    ERR_RETURN
}

void FragLibraryTronDIA::buildChargeVsIonLabels() {

    for (int chrg = m_params.chargeStateMin; chrg <= m_params.chargeStateMax; ++chrg) {
        const QStringList ionLabels = TandemPredictionUtils::buildIonLabels(chrg);
        m_chargeVsIonLabels.insert(chrg, ionLabels);
    }

}

namespace {

    Err buildIonLabelVsDoubleVec(
            const IonLabels &ionLabels,
            const QVector<double> &dVals,
            QHash<IonLabel, double> *ionLabelsVsDVals
    ) {

        ERR_INIT

        e = ErrorUtils::isEqual(dVals.size(),ionLabels.size()); ree;

        for (int i = 0; i < dVals.size(); i++) {

            ionLabelsVsDVals->insert(
                    ionLabels.at(i),
                    dVals.at(i)
            );
        }

        ERR_RETURN
    }

    QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>> buildMS2IonsForPeptideSequenceChargeKey(
            const TandemLibraryReaderRow &tandemLibraryReaderRow,
            const QMap<Charge, IonLabels> &chargeVsIonLabels,
            const AminoAcids &aminoAcids
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(chargeVsIonLabels); rree;

        PeptideString peptideString;
        int charge = -1;
        e = TandemPredictionUtils::extractSequenceAndChargeFromPeptideSequenceChargeKey(
                tandemLibraryReaderRow.peptideSequenceChargeKey,
                &peptideString,
                &charge
        ); rree;

        //TODO write code to extract mods w/i string i.e., PEPT[15.99]IDE
        // For now just use {} in place for modifications.
        QVector<double> mzVals;
        e = TandemPredictionUtils::calculateMzValuesForPrediction(
                peptideString,
                {},
                aminoAcids,
                charge,
                &mzVals
        ); rree;

        const IonLabels ionLabels =  chargeVsIonLabels.value(charge);

        QHash<IonLabel, double> ionLabelsVsMzVals;
        e = buildIonLabelVsDoubleVec(
                ionLabels,
                mzVals,
                &ionLabelsVsMzVals
        ); rree;

        QHash<IonLabel, double> ionLabelsVsIntensityVals;
        e = buildIonLabelVsDoubleVec(
                tandemLibraryReaderRow.ionLabels,
                tandemLibraryReaderRow.intensityVals,
                &ionLabelsVsIntensityVals
        ); rree;

        QVector<MS2Ion> ms2Ions;
        for (auto it = ionLabelsVsIntensityVals.begin(); it != ionLabelsVsIntensityVals.end(); it++) {

            MS2Ion ms2Ion;
            ms2Ion.ionLabel = it.key();
            ms2Ion.intensity = it.value();

            e = ErrorUtils::isTrue(ionLabelsVsMzVals.contains(ms2Ion.ionLabel)); rree;
            ms2Ion.mz = ionLabelsVsMzVals.value(ms2Ion.ionLabel);

            ms2Ions.push_back(ms2Ion);
        }

        const auto sortMzAscLogic = [](const MS2Ion &l, const MS2Ion &r){return l.mz < r.mz;};
        std::sort(ms2Ions.begin(), ms2Ions.end(), sortMzAscLogic);

        return {e , {tandemLibraryReaderRow.peptideSequenceChargeKey, ms2Ions}};
    }


}//namespace
Err FragLibraryTronDIA::readFragLibFile(const QString &fragLibFilePath) {

    ERR_INIT

    QVector<TandemLibraryReaderRow> tandemLibraryRows;
    e = TandemLibraryReader::readTandemPredictions(
            fragLibFilePath,
            &tandemLibraryRows
    );

    qDebug() << "Library size:" << tandemLibraryRows.size();

#define RUN_PARALLEL_MZ_ASS
#ifdef RUN_PARALLEL_MZ_ASS

    qDebug() << "Running FragLibraryTronDIA parallel";
        const auto mzBuilderLogicBinder = std::bind(
                buildMS2IonsForPeptideSequenceChargeKey,
                std::placeholders::_1,
                m_chargeVsIonLabels,
                m_params.aminoAcids
        );

        QFuture<QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>>> futures = QtConcurrent::mapped(
                tandemLibraryRows,
                mzBuilderLogicBinder
        );
        futures.waitForFinished();

        for (const QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>> &pepMS2Ions : futures) {

            e = pepMS2Ions.first; ree;
            const QPair<PeptideSequenceChargeKey , QVector<MS2Ion>> &res = pepMS2Ions.second;
            m_pepSeqChrgKeyVsMS2Ions.insert(res.first, res.second);
        }
#else
    qDebug() << "Running FragLibraryTronDIA serial";
    for (const TandemLibraryReaderRow &tandemLibraryReaderRow : tandemLibraryRows) {

        QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>> pepMS2Ions
                = buildMS2IonsForPeptideSequenceChargeKey(
                        tandemLibraryReaderRow,
                        m_chargeVsIonLabels,
                        m_params.aminoAcids
                        );
        e = pepMS2Ions.first; ree;

        const QPair<PeptideSequenceChargeKey , QVector<MS2Ion>> &res = pepMS2Ions.second;
        m_pepSeqChrgKeyVsMS2Ions.insert(res.first, res.second);
    }
#endif

    ERR_RETURN
}




