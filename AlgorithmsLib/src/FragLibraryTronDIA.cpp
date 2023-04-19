//
// Created by anichols on 4/10/23.
//

#include "FragLibraryTronDIA.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "TandemFragmentPredictotron.h"
#include "TandemLibraryReader.h"
#include "TandemPredictionUtils.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

#include <iostream>

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
    qDebug() << "Peptides count charge vs uniques" << m_pepSeqChrgKeyVsMS2Ions.size() << m_peptideWithModsVsisDecoy.size();

    ERR_RETURN
}

namespace {

    const auto sortMzAscLogic = [](const MS2Ion &l, const MS2Ion &r){return l.mz < r.mz;};

    const auto sortIntzAscLogic
            = [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;};

    QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>> buildMS2IonsForPeptideSequenceChargeKey(
            const TandemLibraryReaderRow &tandemLibraryReaderRow
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(tandemLibraryReaderRow.mzVals); rree;
        e = ErrorUtils::isEqual(
                tandemLibraryReaderRow.mzVals.size(),
                tandemLibraryReaderRow.intensityVals.size()
                ); rree;

        e = ErrorUtils::isEqual(
                tandemLibraryReaderRow.mzVals.size(),
                tandemLibraryReaderRow.ionLabels.size()
        ); rree;

        PeptideString peptideString;
        int charge = -1;
        e = TandemPredictionUtils::extractSequenceAndChargeFromPeptideSequenceChargeKey(
                tandemLibraryReaderRow.peptideSequenceChargeKey,
                &peptideString,
                &charge
        ); rree;

        QVector<MS2Ion> ms2Ions;
        for (int i = 0; i < tandemLibraryReaderRow.mzVals.size(); i++) {

            MS2Ion ms2Ion;
            ms2Ion.ionLabel = tandemLibraryReaderRow.ionLabels.at(i);
            ms2Ion.intensity = tandemLibraryReaderRow.intensityVals.at(i);
            ms2Ion.mz = tandemLibraryReaderRow.mzVals.at(i);

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
    ); ree;

    qDebug() << "Library size:" << tandemLibraryRows.size();

#define RUN_PARALLEL_MZ_ASS
#ifdef RUN_PARALLEL_MZ_ASS
    qDebug() << "Running FragLibraryTronDIA parallel";
    QFuture<QPair<Err, QPair<PeptideSequenceChargeKey , QVector<MS2Ion>>>> futures = QtConcurrent::mapped(
            tandemLibraryRows,
            buildMS2IonsForPeptideSequenceChargeKey
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
                = buildMS2IonsForPeptideSequenceChargeKey(tandemLibraryReaderRow);
        e = pepMS2Ions.first; ree;

        const QPair<PeptideSequenceChargeKey , QVector<MS2Ion>> &res = pepMS2Ions.second;
        m_pepSeqChrgKeyVsMS2Ions.insert(res.first, res.second);
    }
#endif

    e = buildPeptideWithModsVsisDecoy(tandemLibraryRows); ree;

    ERR_RETURN
}

Err FragLibraryTronDIA::buildPeptideWithModsVsisDecoy(
        const QVector<TandemLibraryReaderRow> &tandemLibraryRows
        ) {

    ERR_INIT

    for (const TandemLibraryReaderRow &row : tandemLibraryRows) {

        const PeptideSequenceChargeKey &peptideSequenceChargeKey = row.peptideSequenceChargeKey;
        const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                QString::SkipEmptyParts
                );

        const int expectedSplitSize = 2;
        e = ErrorUtils::isEqual(
                peptideSequenceChargeKeySplit.size(),
                expectedSplitSize
                ); ree;

        const PeptideStringWithMods &peptideStringWithMods = peptideSequenceChargeKeySplit.front();

        m_peptideWithModsVsisDecoy.insert(
                peptideStringWithMods,
                row.isDecoy
                );
    }

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

    ms2Ions->clear();

    e = ErrorUtils::isTrue(m_pepSeqChrgKeyVsMS2Ions.contains(peptideSequenceChargeKey));
    if (e != eNoError) {
        const QString pep = peptideSequenceChargeKey.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP).at(0);
        qDebug() << "Prediction not found" << peptideSequenceChargeKey
                 << BiophysicalCalcs::calculatePeptideMass(pep, m_params.aminoAcids, {});
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

    const double mzMin = 0.0;
    const double mzMax = 2000.0;

    e = getMS2Ions(
            peptideSequenceChargeKey,
            {mzMin, mzMax},
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

    peptideStringWithModsVsMS2Ions->clear();
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
            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
                    peptideStringWithMods,
                    charge
            );


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

bool FragLibraryTronDIA::isInit() {
    return !m_pepSeqChrgKeyVsMS2Ions.isEmpty();
}

Err FragLibraryTronDIA::peptideStringWithModsIsDecoy(
        const PeptideStringWithMods &peptideStringWithMods,
        bool *isDecoy
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptideWithModsVsisDecoy); ree;
    e = ErrorUtils::isTrue(m_peptideWithModsVsisDecoy.contains(peptideStringWithMods)); ree;
    *isDecoy = m_peptideWithModsVsisDecoy.value(peptideStringWithMods);

    ERR_RETURN
}

ScanPoints FragLibraryTronDIA::ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions) {

    ScanPoints scanPoints;
    for (const MS2Ion &ion : ms2Ions) {
        scanPoints.push_back({ion.mz, ion.intensity});
    }

    return scanPoints;
}
