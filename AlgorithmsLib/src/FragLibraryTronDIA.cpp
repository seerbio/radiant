//
// Created by anichols on 4/10/23.
//

#include "FragLibraryTronDIA.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "TandemLibraryReader.h"
#include "TandemPredictionUtils.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

namespace {

    void stripToPeptideStringWithMods(
            const QList<PeptideSequenceChargeKey> &peptideSequenceChargeKey,
            QVector<PeptideStringWithMods> *peptideStringWithMods
            ) {
        const auto transformLogic = [&](const PeptideSequenceChargeKey &pk) {
            return pk.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP).at(0);
        };

        std::transform(
                peptideSequenceChargeKey.begin(),
                peptideSequenceChargeKey.end(),
                std::back_inserter(*peptideStringWithMods),
                transformLogic
        );
    }

}//namespace
Err FragLibraryTronDIA::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibFilePath
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    m_params = pythiaParameters;

    e = buildChargeVsIonLabels(); ree;

    e = readFragLibFile(fragLibFilePath); ree;

    const QList<PeptideSequenceChargeKey> &peptideSequenceChargeKey = m_pepSeqChrgKeyVsMS2Ions.keys();
    QVector<PeptideStringWithMods> peptideStringWithMods;
    stripToPeptideStringWithMods(
            peptideSequenceChargeKey,
            &peptideStringWithMods
            );

    e = m_peptideMassRTree.init(
            peptideStringWithMods,
            m_params.aminoAcids
            ); ree;

    qDebug() << "Library init in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err FragLibraryTronDIA::buildChargeVsIonLabels() {

    ERR_INIT

    for (int chrg = m_params.chargeStateMin; chrg <= m_params.chargeStateMax; ++chrg) {
        const QStringList ionLabels = TandemPredictionUtils::buildIonLabels(chrg);
        m_chargeVsIonLabels.insert(chrg, ionLabels);
    }

    e = ErrorUtils::isNotEmpty(m_chargeVsIonLabels); ree;

    ERR_RETURN
}

namespace {

    const auto sortMzAscLogic = [](const MS2Ion &l, const MS2Ion &r){return l.mz < r.mz;};

    const auto sortIntzAscLogic
            = [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;};

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

        std::sort(ms2Ions.begin(), ms2Ions.end(), sortMzAscLogic);

        return {e , {tandemLibraryReaderRow.peptideSequenceChargeKey, ms2Ions}};
    }


}//namespace
Err FragLibraryTronDIA::readFragLibFile(const QString &fragLibFilePath) {

    ERR_INIT

    QVector<TandemLibraryReaderRow> tandemLibraryRows;
    e = ParquetReader::read(
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

namespace {

    void filterMS2IonsByMz(
            const QPair<double, double> mzStartMzEnd,
            QVector<MS2Ion> *ms2Ions
            ) {

        const auto terminatorLogic = [mzStartMzEnd](const MS2Ion &ms2Ion){
            return ms2Ion.mz < mzStartMzEnd.first || ms2Ion.mz > mzStartMzEnd.second;
        };

        const auto terminator = std::remove_if(ms2Ions->begin(), ms2Ions->end(), terminatorLogic);
        ms2Ions->erase(terminator, ms2Ions->end());
    }

    void filterMs2IonsTopNIntense(
            int topN,
            QVector<MS2Ion> *ms2Ions
            ) {

        topN = std::min(topN, ms2Ions->size());

        std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntzAscLogic);

        ms2Ions->resize(topN);

        std::sort(ms2Ions->begin(), ms2Ions->end(), sortMzAscLogic);
    }

}//namespace
Err FragLibraryTronDIA::getMS2Ions(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        QPair<double, double> mzStartMzEnd,
        QVector<MS2Ion> *ms2Ions
        ) {

    ERR_INIT

    if (mzStartMzEnd.first > 989) {
        qDebug() << peptideSequenceChargeKey;
    }

    e = ErrorUtils::isTrue(m_pepSeqChrgKeyVsMS2Ions.contains(peptideSequenceChargeKey));
    if (e != eNoError) {
        const QString pep = peptideSequenceChargeKey.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP).at(0);
        qDebug() << peptideSequenceChargeKey << BiophysicalCalcs::calculatePeptideMass(pep, m_params.aminoAcids, {});
        rrr(eValueError);
    }

    *ms2Ions = m_pepSeqChrgKeyVsMS2Ions.value(peptideSequenceChargeKey);

    filterMS2IonsByMz(
            mzStartMzEnd,
            ms2Ions
            );

    ERR_RETURN
}

Err FragLibraryTronDIA::getMS2Ions(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        QVector<MS2Ion> *ms2Ions
        ) {

    ERR_INIT
    e = getMS2Ions(
            peptideSequenceChargeKey,
            {0.0, 2000.0},
            ms2Ions
            ); ree;
    ERR_RETURN
}

Err FragLibraryTronDIA::getMS2Ions(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        QPair<double, double> mzStartMzEnd,
        int topNIntense,
        QVector<MS2Ion> *ms2Ions
        ) {

    ERR_INIT

    e = getMS2Ions(
            peptideSequenceChargeKey,
            mzStartMzEnd,
            ms2Ions
    ); ree;

    if (topNIntense > 0) {
        filterMs2IonsTopNIntense(
                topNIntense,
                ms2Ions
        );
    }

    ERR_RETURN
}

Err FragLibraryTronDIA::getMS2Ions(
        double mzTargetStart,
        double mzTargetEnd,
        int topNIntense,
        QMap<PeptideStringWithMods, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions
        ) {

    ERR_INIT

    const int monoOffset = 0;

    for (int charge = m_params.chargeStateMin; charge <= m_params.chargeStateMax; charge++) {

        const double massStart = BiophysicalCalcs::calculateMassFromThomson(
                mzTargetStart,
                charge,
                monoOffset
                );

        const double massEnd = BiophysicalCalcs::calculateMassFromThomson(
                mzTargetEnd,
                charge,
                monoOffset
        );

        QHash<PeptideStringWithMods , Mass> peptideStringWithModsTableVsMass;
        e = m_peptideMassRTree.getPeptides(
                massStart,
                massEnd,
                &peptideStringWithModsTableVsMass
        ); ree;

        for (auto it = peptideStringWithModsTableVsMass.begin();
                  it != peptideStringWithModsTableVsMass.end();
                  it++
                  ) {

            const Mass mass = it.value();
            const double mzTargetTheo = BiophysicalCalcs::calculateThomsonFromMass(mass, charge);
            const PeptideStringWithMods  &peptideStringWithMods = it.key();
            const PeptideSequenceChargeKey peptideSequenceChargeKey
                = peptideStringWithMods + S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP + QString::number(charge);

            if (mzTargetTheo < mzTargetStart || mzTargetTheo > mzTargetEnd) {
                continue;
            }

            QVector<MS2Ion> ms2Ions;
            e = getMS2Ions(
                    peptideSequenceChargeKey,
                    {m_params.mzMinDataStructure, m_params.mzMaxDataStructure},
                    topNIntense,
                    &ms2Ions
                    ); ree;

            peptideStringWithModsVsMS2Ions->insert(peptideStringWithMods, ms2Ions);
        }
    }

    ERR_RETURN
}
