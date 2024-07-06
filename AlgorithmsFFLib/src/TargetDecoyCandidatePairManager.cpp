//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePairManager.h"

#include "AminoAcids.h"
#include "BiophysicalCalcs.h"
#include "FragLibReader.h"
#include "MolecularFormula.h"
#include "Molecule.h"
#include "ParallelUtils.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>


Err TargetDecoyCandidatePairManager::init(
        const PythiaParameters &pythiaParameters,
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isFalse(fragLibReaderRows->isEmpty()); ree;

    m_pythiaParameters = pythiaParameters;

    e = buildTargetDecoyCandidatePairs(*fragLibReaderRows); ree;

    ERR_RETURN
}

namespace {

    Err buildMS2Ions(
            const PythiaParameters &pythiaParameters,
            const FragLibReaderRow &flrr,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(flrr.mzVals); ree;
        e = ErrorUtils::isEqual(flrr.mzVals.size(), flrr.intensityVals.size());ree;
        e = ErrorUtils::isNotEmpty(flrr.ionLabels); ree;

        QString ionLabels = flrr.ionLabels;

        const QStringList ionLabelsSplit = ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, Qt::SkipEmptyParts);
        e = ErrorUtils::isEqual(flrr.mzVals.size(), ionLabelsSplit.size());ree;

        QVector<MS2Ion> ms2IonsBuilder;
        ms2IonsBuilder.reserve(flrr.mzVals.size());
        for (int i = 0; i < flrr.mzVals.size(); i++) {
            MS2Ion ms2Ion;
            ms2Ion.mz = flrr.mzVals.at(i);
            ms2Ion.intensity = flrr.intensityVals.at(i);
            ms2Ion.ionLabel = ionLabelsSplit.at(i);
            ms2Ion.charge = ms2Ion.ionLabel.contains("^2") ? 2 : 1;

            ms2IonsBuilder.push_back(ms2Ion);
        }

        MS2Ion::filterMS2IonsByMz(
                pythiaParameters.mzMinMS2,
                pythiaParameters.mzMaxMS2,
                &ms2IonsBuilder
                );

        MS2Ion::sortMS2IonsIntensityDesc(&ms2IonsBuilder);
        for (int i = 0; i < ms2IonsBuilder.size(); i++) {
            MS2Ion &ms2Ion = ms2IonsBuilder[i];
            ms2Ion.rank = i;
        }

        const int ms2IonsSize = ms2IonsBuilder.size();
        ms2Ions->reserve(ms2IonsSize);
        ms2IonsBuilder.resize(ms2IonsSize);

        *ms2Ions = ms2IonsBuilder;

        ERR_RETURN
    }

    Err mutateCandidatePeptideTarget(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonTarget,
            QVector<MS2Ion> *ms2IonDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;
        e = ErrorUtils::isNotEmpty(ms2IonTarget); ree;

        ms2IonDecoys->clear();
        ms2IonDecoys->reserve(ms2IonTarget.size());

        const QMap<QChar, double> diannMutateAminoAcidToMass = AminoAcids::diannMutateAminoAcidToMass();

        const PeptideString &peptideString = peptideStringWithMods.removeUniModChars();

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = peptideString.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass; //NOTE: do not static cast to float
                    }
                }
            }

            else if (ionLableInfo.second.contains('y')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass; //NOTE: do not static cast to float
                    }
                }
            }

            else {
                qDebug() << "Non b/y/a ion" << ionLableInfo;
            }

            ms2IonDecoys->push_back(ms2IonDecoy);
        }

        ERR_RETURN
    }

    Err reverseCandidatePeptideTarget(
            const PeptideStringWithMods &peptideStringWithMods,
            const AminoAcids &aminoAcids,
            const QVector<MS2Ion> &ms2IonTarget,
            QVector<MS2Ion> *ms2IonDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;
        e = ErrorUtils::isNotEmpty(ms2IonTarget); ree;

        ms2IonDecoys->clear();
        ms2IonDecoys->reserve(ms2IonTarget.size());

        PeptideStringWithMods peptideStringWithModsMiddelReversed = peptideStringWithMods;
        std::reverse(
                peptideStringWithModsMiddelReversed.begin() + 1,
                peptideStringWithModsMiddelReversed.begin() + peptideStringWithModsMiddelReversed.size() - 1
                );

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b')) {
                const QString bSeq = peptideStringWithModsMiddelReversed.left(ionLableInfo.first);
                ms2IonDecoy.mz = static_cast<float>(BiophysicalCalcs::calculateThomson(bSeq, aminoAcids, ms2IonDecoy.charge));
            }

            else if (ionLableInfo.second.contains('y')) {

                const QString ySeq = peptideStringWithModsMiddelReversed.right(ionLableInfo.first);
                ms2IonDecoy.mz
                    = static_cast<float>(BiophysicalCalcs::calculateThomson(ySeq, aminoAcids, ms2IonDecoy.charge) +  (CommonMolecules::H2O.monoisotopicMass() / ms2Ion.charge));
            }

            else {
                qDebug() << "Non b/y/a ion" << ionLableInfo;
            }

            ms2IonDecoys->push_back(ms2IonDecoy);
        }

        ERR_RETURN
    }


    QPair<Err, TargetDecoyCandidatePair> targetDecoyCandidatePairsLoadLogic(
            const FragLibReaderRow &flrr,
            const PythiaParameters &pythiaParameters
            ) {

        ERR_INIT

        PeptideStringWithMods peptideStringWithMods;
        Charge charge;
        e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
                flrr.peptideSequenceChargeKey,
                &peptideStringWithMods,
                &charge
        ); rree;

        const int peptideLength = peptideStringWithMods.sizeNoMods();
        if (peptideLength < pythiaParameters.peptideLengthMin
            || peptideLength > pythiaParameters.peptideLengthMax
            || charge < pythiaParameters.chargeStateMin
            || charge > pythiaParameters.chargeStateMax
            || flrr.isDecoy == 1) {
            return {e, {}};
        }

        QVector<MS2Ion> ms2IonsTarget;
        e = buildMS2Ions(
                pythiaParameters,
                flrr,
                &ms2IonsTarget
        ); rree;

        if (ms2IonsTarget.isEmpty()) {
            return{eValueError, {}};
        }

        QVector<MS2Ion> ms2IonsDecoy;

        e= mutateCandidatePeptideTarget(
                peptideStringWithMods,
                ms2IonsTarget,
                &ms2IonsDecoy
        ); rree;

        TargetDecoyCandidatePair targetDecoyCandidatePair(
                peptideStringWithMods,
                ms2IonsTarget,
                ms2IonsDecoy,
                flrr.precursorCharge,
                static_cast<float>(flrr.mass),
                static_cast<float>(flrr.iRT),
                flrr.mzVals.size()
        );

        return {e, targetDecoyCandidatePair};
    }

}//namespace
Err TargetDecoyCandidatePairManager::buildTargetDecoyCandidatePairs(
        const QVector<FragLibReaderRow> &fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibReaderRows); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;

    QElapsedTimer et;
    et.start();

    m_targetDecoyCandidatePairs.reserve(fragLibReaderRows.size()); ree;

#define PARALLEL_FRAGLIB_LOAD
#ifdef PARALLEL_FRAGLIB_LOAD
    const auto loadLogicBinder = std::bind(
            targetDecoyCandidatePairsLoadLogic,
            std::placeholders::_1,
            m_pythiaParameters
    );

    QFuture<QPair<Err, TargetDecoyCandidatePair>> futures = QtConcurrent::mapped(
            fragLibReaderRows,
            loadLogicBinder
            );
    futures.waitForFinished();

    for (const QPair<Err, TargetDecoyCandidatePair> &result : futures) {
        if (result.first == eValueError) {
            continue;
        }
        e = result.first; ree;
        const TargetDecoyCandidatePair &tdcp = result.second;
        m_targetDecoyCandidatePairs.push_back(tdcp);
    }
#else
    for (const FragLibReaderRow &flrr : fragLibReaderRows) {

        QPair<Err, TargetDecoyCandidatePair> result = targetDecoyCandidatePairsLoadLogic(
                flrr,
                m_pythiaParameters
                );

        if (result.first == eValueError) {
            continue;
        }

        e = result.first; ree;

        const TargetDecoyCandidatePair &tdcp = result.second;
        m_targetDecoyCandidatePairs.push_back(tdcp);
    }
#endif

    e = filterDecoySequencesThatAreAlsoTargetSequences();

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << m_targetDecoyCandidatePairs.size() << "Candidates loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    int mangleLogic(
        const QHash<PeptideString, bool> &peptideStringIsoleucineReplaceVsIsAlsoDecoy,
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairsPntrs
        ) {

        int modified = 0;
        for (TargetDecoyCandidatePair *tdcp : targetDecoyCandidatePairsPntrs) {

            const PeptideString peptideStringWithModsMutated
                = AminoAcids::mutatePenultimatePeptideResidues(tdcp->peptideStringWithMods()).replace('I', 'L');

            if (peptideStringIsoleucineReplaceVsIsAlsoDecoy.contains(peptideStringWithModsMutated)) {
                tdcp->mangleMs2IonsDecoy();
                modified++;
            }
        }

        return modified;
    }

}
Err TargetDecoyCandidatePairManager::filterDecoySequencesThatAreAlsoTargetSequences() {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;

    QElapsedTimer et;
    et.start();

    QHash<PeptideString, bool> peptideStringIsoleucineReplaceVsIsAlsoDecoy;
    for (TargetDecoyCandidatePair &tdcp : m_targetDecoyCandidatePairs) {
        peptideStringIsoleucineReplaceVsIsAlsoDecoy.insert(tdcp.peptideStringWithMods().replace('I', 'L'), false);
    }

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs;
    e = getTargetDecoyCandidatePairPointers(&targetDecoyCandidatePairsPntrs); ree;

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        targetDecoyCandidatePairsPntrs,
        m_pythiaParameters.threadCount,
        &targetDecoyCandidatePairsPntrsTranched
        ); ree;

    const auto manglerLogicBinder = std::bind(
        mangleLogic,
        peptideStringIsoleucineReplaceVsIsAlsoDecoy,
        std::placeholders::_1
        );

    QFuture<int> futures = QtConcurrent::mapped(
        targetDecoyCandidatePairsPntrsTranched,
        manglerLogicBinder
        );
    futures.waitForFinished();

    int modified = 0;
    for (int result : futures) {
        modified += result;
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << modified << "Sequences were found to have decoys that were also targets and were modified!!!!" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairsPntrs) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;

    for (TargetDecoyCandidatePair &t : m_targetDecoyCandidatePairs) {

        const int peptideStringsize = t.peptideString().size();
        if (!(m_pythiaParameters.peptideLengthMin <= peptideStringsize && peptideStringsize <= m_pythiaParameters.peptideLengthMax)) {
            continue;
        }

        targetDecoyCandidatePairsPntrs->push_back(&t);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        PeptideStringWithMods *peptideStringWithMods,
        Charge *charge
        ){

    ERR_INIT

    const int expectedSplitSize = 2;

    const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
            S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
            Qt::SkipEmptyParts
    );

    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize
            ); ree;

    *peptideStringWithMods = PeptideStringWithMods(peptideSequenceChargeKeySplit.front());

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
    ); ree

    ERR_RETURN
}

bool TargetDecoyCandidatePairManager::isInit() const {
    return !m_targetDecoyCandidatePairs.empty();
}

int TargetDecoyCandidatePairManager::targetsCount() const {
    return m_targetDecoyCandidatePairs.size();
}
