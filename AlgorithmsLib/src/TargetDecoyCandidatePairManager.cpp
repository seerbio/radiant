//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePairManager.h"

#include "AminoAcids.h"
#include "FragLibReader.h"
#include "MolecularFormula.h"
#include "Molecule.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>


TargetDecoyCandidatePairManager::TargetDecoyCandidatePairManager() : m_isInit(false) {}

bool TargetDecoyCandidatePairManager::isInit() {
    return m_isInit;
}


Err TargetDecoyCandidatePairManager::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibFileUri
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibFileUri); ree;

    m_pythiaParameters = pythiaParameters;

    const double massMin
        = pythiaParameters.peptideLengthMin * Molecule(MolecularFormulas::alanineFormula).monoisotopicMass();

    const double massMax
            = pythiaParameters.peptideLengthMax * Molecule(MolecularFormulas::tryptophanFormula).monoisotopicMass();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            fragLibFileUri,
            massMin,
            massMax,
            &fragLibReaderRows
            ); ree;

    e = buildTargetDecoyCandidatePairs(fragLibReaderRows); ree;
    e = buildIndexVsTargetDecoyCandidatePairPtrs(); ree;
    e = initBoostRTreeWrapper(); ree;

    m_isInit = true;

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
        e = ErrorUtils::isNotEmpty(flrr.ionLabels); ree

        const QStringList ionLabelsSplit = flrr.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);
        e = ErrorUtils::isEqual(flrr.mzVals.size(), ionLabelsSplit.size());ree;

        QVector<MS2Ion> ms2IonsBuilder;
        ms2IonsBuilder.reserve(flrr.mzVals.size());

        for (int i = 0; i < flrr.mzVals.size(); i++) {
            MS2Ion ms2Ion;
            ms2Ion.mz = flrr.mzVals.at(i);
            ms2Ion.intensity = flrr.intensityVals.at(i);
            ms2Ion.ionLabel = ionLabelsSplit.at(i);
            ms2Ion.iRT = static_cast<IRT>(flrr.iRT);
            ms2Ion.charge = flrr.charge;

            ms2IonsBuilder.push_back(ms2Ion);
        }

        MS2Ion::filterMS2IonsByMz(
                pythiaParameters.mzMinDataStructure,
                pythiaParameters.mzMaxDataStructure,
                &ms2IonsBuilder
        );

        MS2Ion::sortMS2IonsIntensityDesc(&ms2IonsBuilder);

        const int ms2IonsSize = std::min(pythiaParameters.topNMs2Ions, ms2IonsBuilder.size());
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

        const QMap<QChar, double> diannMutateAminoAcidTo = AminoAcids::diannMutateAminoAcidTo();

        const QString &seq = peptideStringWithMods;

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = seq.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }
                }
            }

            else if (ionLableInfo.second.contains('y')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
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

        //TODO once incorporating PTMs in the sequence, make algo to remove them
        // and then calculate the peptide length MODS
        const int peptideLength = peptideStringWithMods.size();
        if (peptideLength < pythiaParameters.peptideLengthMin
            || peptideLength > pythiaParameters.peptideLengthMax
            || flrr.isDecoy) {
            return {e, {}};
        }

        QVector<MS2Ion> ms2IonsTarget;
        e = buildMS2Ions(
                pythiaParameters,
                flrr,
                &ms2IonsTarget
        ); rree;

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
                flrr.charge,
                flrr.mass,
                flrr.iRT,
                flrr.mzVals.size(),
                -1
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

        e = result.first; ree;

        const TargetDecoyCandidatePair &tdcp = result.second;
        m_targetDecoyCandidatePairs.push_back(tdcp);
    }
#endif

    qDebug() << m_targetDecoyCandidatePairs.size() << "Candidates loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::buildIndexVsTargetDecoyCandidatePairPtrs() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;

    std::sort(
            m_targetDecoyCandidatePairs.begin(),
            m_targetDecoyCandidatePairs.end(),
            [](const TargetDecoyCandidatePair &l, const TargetDecoyCandidatePair &r){return l.mz() < r.mz();}
            );

    for (TargetDecoyCandidatePair &tdcp : m_targetDecoyCandidatePairs) {

        const int currentIndexVsTargetDecoyCandidatePairPtrsSize = m_indexVsTargetDecoyCandidatePairPtrs.size();

        tdcp.setTargetDecoyCandidatePairIndex(currentIndexVsTargetDecoyCandidatePairPtrsSize);
        m_indexVsTargetDecoyCandidatePairPtrs.insert(
                currentIndexVsTargetDecoyCandidatePairPtrsSize,
                &tdcp
                );
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::initBoostRTreeWrapper() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isNotEmpty(m_indexVsTargetDecoyCandidatePairPtrs); ree;

    QElapsedTimer et;
    et.start();

    QVector<RTreePointData2D> dataPoints;
    for (const TargetDecoyCandidatePair &tdcp : m_targetDecoyCandidatePairs) {
        dataPoints.push_back({tdcp.mz(), 0.0, static_cast<double>(tdcp.targetDecoyCandidatePairIndex())});
    }

    e = m_boostRTreeWrapper.init(dataPoints); ree;

    qDebug() <<"TargetDecoyPairs RTree loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(
        double mzMin,
        double mzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isNotEmpty(m_indexVsTargetDecoyCandidatePairPtrs); ree;
    e = ErrorUtils::isTrue(mzMax >= mzMin); ree;

    QVector<RTreePointData2D> vals;
    e = m_boostRTreeWrapper.getPoints(
            mzMin,
            mzMax,
            0.0,
            0.0,
            &vals
            ); ree;

    const auto transformLogic = [&](const RTreePointData2D &p){
        return m_indexVsTargetDecoyCandidatePairPtrs.value(static_cast<TargetDecoyCandidatePairIndex>(p.val));
    };

    std::transform(
            vals.begin(),
            vals.end(),
            std::back_inserter(*targetDecoyPointers),
            transformLogic
            );

    std::mt19937 gen(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
    std::shuffle(targetDecoyPointers->begin(), targetDecoyPointers->end(), gen);

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(
        double mzMin,
        double mzMax,
        double randomSelectionFraction,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

    ERR_INIT

    QVector<TargetDecoyCandidatePair*> targetDecoyPointersAll;
    e = getTargetDecoyCandidatePairPointers(
            mzMin,
            mzMax,
            &targetDecoyPointersAll
            ); ree;

    const int testDataSize = static_cast<int>(targetDecoyPointersAll.size() * randomSelectionFraction);
    const int seed = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;

    const QMap<int, bool> selectionList = MathUtils::generateRandomSelectionList(
            targetDecoyPointersAll.size(),
            testDataSize,
            seed
            );

    for (int i = 0; i < targetDecoyPointersAll.size(); i++) {

        if (!selectionList.value(i)) {
            continue;
        }

        targetDecoyPointers->push_back(targetDecoyPointersAll[i]);
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
            QString::SkipEmptyParts
    );

    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize); ree;

    *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
    ); ree

    ERR_RETURN
}
