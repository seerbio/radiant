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

    e = buildTargetDecoyCandidatePairs(fragLibReaderRows); ree;

    ERR_RETURN
}

namespace {

    float calculateDecoyMassDelta(const PeptideStringWithMods &peptideStringWithMods) {

        const QMap<QChar, double> diannMutateAminoAcidToMass = AminoAcids::diannMutateAminoAcidToMass();

        const PeptideString &peptideString = peptideStringWithMods.removeUniModChars();

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = peptideString.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(secondIndexToMutate));

        return nTermDeltaMass + cTermDeltaMass;
    }

    void targetDecoyCandidatePairsLoadLogic(TargetDecoyCandidatePair &tdcp) {



		tdcp.init();

    }

    void filterNullSequences(QVector<TargetDecoyCandidatePair> *targetDecoyCandidatePairs) {
        const auto terminatorLogic = [](const TargetDecoyCandidatePair &tdcp){return tdcp.peptideStringWithMods().isEmpty();};
        const auto terminator = std::remove_if(targetDecoyCandidatePairs->begin(), targetDecoyCandidatePairs->end(), terminatorLogic);
        targetDecoyCandidatePairs->erase(terminator, targetDecoyCandidatePairs->end());
    }

    void filterLeucineIsoleucineIsomers(QVector<TargetDecoyCandidatePair> *targetDecoyCandidatePairs) {

        QVector<TargetDecoyCandidatePair> targetDecoyCandidatePairsFiltered;

        QHash<QString, bool> isomers;
        for (int i = 0; i < targetDecoyCandidatePairs->size(); ++i) {

            const TargetDecoyCandidatePair &tdcp = targetDecoyCandidatePairs->at(i);

            QString peptideStringWithModsIsomer = tdcp.peptideStringWithMods();
            peptideStringWithModsIsomer = peptideStringWithModsIsomer.replace('I', 'L');
            peptideStringWithModsIsomer += QString::number(tdcp.charge());

            if (isomers.value(peptideStringWithModsIsomer)) {
                continue;
            }

            targetDecoyCandidatePairsFiltered.push_back(tdcp);
            isomers.insert(peptideStringWithModsIsomer, true);
        }

        *targetDecoyCandidatePairs = targetDecoyCandidatePairsFiltered;

    }

}//namespace
Err TargetDecoyCandidatePairManager::buildTargetDecoyCandidatePairs(
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(fragLibReaderRows->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;

    QElapsedTimer et;
    et.start();

    m_targetDecoyCandidatePairs.resize(fragLibReaderRows->size()); ree;
	e = ErrorUtils::isEqual(fragLibReaderRows->size(), m_targetDecoyCandidatePairs.size()); ree;
	for (int i = 0; i < m_targetDecoyCandidatePairs.size(); i++) {
		m_targetDecoyCandidatePairs[i].setFragLibReaderRowPntr(&(*fragLibReaderRows)[i]);
	}

#define PARALLEL_FRAGLIB_LOAD
#ifdef PARALLEL_FRAGLIB_LOAD

    QFuture<void> futures = QtConcurrent::map(
            m_targetDecoyCandidatePairs,
            targetDecoyCandidatePairsLoadLogic
            );
    futures.waitForFinished();


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

    filterNullSequences(&m_targetDecoyCandidatePairs);
    filterLeucineIsoleucineIsomers(&m_targetDecoyCandidatePairs);

    e = filterDecoySequencesThatAreAlsoTargetSequences(); ree;

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
                = AminoAcids::mutatePeptideResidues(tdcp->peptideStringWithMods(), 1).replace('I', 'L');

            if (peptideStringIsoleucineReplaceVsIsAlsoDecoy.contains(peptideStringWithModsMutated)) {
                tdcp->decoySharesSequenceWithOtherTarget(true);
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

namespace {

    void logic(
        PythiaParameters pythiaParameters,
        QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs
        ) {

        const auto terminatorLogic = [pythiaParameters](TargetDecoyCandidatePair *tdcp) {
            const int peptideStringsize = tdcp->peptideString().size();
            return !(pythiaParameters.peptideLengthMin <= peptideStringsize && peptideStringsize <= pythiaParameters.peptideLengthMax);
        };
        const auto terminator = std::remove_if(
            targetDecoyCandidatePairs.begin(),
            targetDecoyCandidatePairs.end(),
            terminatorLogic
            );

        targetDecoyCandidatePairs.erase(terminator, targetDecoyCandidatePairs.end());
    }

}//namespace
Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairsPntrs) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;
    for (TargetDecoyCandidatePair &t : m_targetDecoyCandidatePairs) {
        targetDecoyCandidatePairsPntrs->push_back(&t);
    }

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        *targetDecoyCandidatePairsPntrs,
        m_pythiaParameters.threadCount,
        &targetDecoyCandidatePairsPntrsTranched
        ); ree;

    const auto binder = std::bind(
        logic,
        m_pythiaParameters,
        std::placeholders::_1
        );

    QFuture<void> futures = QtConcurrent::map(
        targetDecoyCandidatePairsPntrsTranched,
        binder
        );
    futures.waitForFinished();

    ERR_RETURN
}

bool TargetDecoyCandidatePairManager::isInit() const {
    return !m_targetDecoyCandidatePairs.empty();
}

int TargetDecoyCandidatePairManager::targetsCount() const {
    return m_targetDecoyCandidatePairs.size();
}
