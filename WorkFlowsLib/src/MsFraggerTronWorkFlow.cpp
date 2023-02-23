//
// Created by anichols on 2/19/23.
//

#include "MsFraggerTronWorkFlow.h"

#include "FragLibraryTron.h"
#include "MsReaderMzML.h"
#include "ParallelUtils.h"
#include "PeptidesLibraryTron.h"

#include <QtConcurrent/QtConcurrent>
#include <QFuture>

#include <iostream>


Err MsFraggerTronWorkFlow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &pepLibUri
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibUri); ree;
    e = ErrorUtils::isNotEmpty(pepLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;
    m_pepLibUri = pepLibUri;

    ERR_RETURN
}

namespace {

    void filterTandemScanIonsByMz(
            MzMin mzMin,
            MzMax mzMax,
            QVector<TandemScanIon> *tandemScanIons
            ) {

        const auto terminatorLogic = [mzMin, mzMax](const TandemScanIon &tsi){
            return !(mzMin <= tsi.mz && tsi.mz <= mzMax);
        };

        const auto terminator = std::remove_if(
                tandemScanIons->begin(),
                tandemScanIons->end(),
                terminatorLogic
                );

        tandemScanIons->erase(terminator, tandemScanIons->end());
    }

    Err trancheTandemScanIons(
            const PythiaParameters &pythiaParameters,
            const QString &mzmLFileURI,
            QMap<int, QVector<TandemScanIon>> *tranchedTandemScanIons
            ) {

        ERR_INIT

        MsReaderMzML mzMLReader;
        e = mzMLReader.openFile(mzmLFileURI); ree;

        QVector<TandemScanIon> tandemScanIons = mzMLReader.tandemScanIons();
        filterTandemScanIonsByMz(
                pythiaParameters.mzMinDataStructure,
                pythiaParameters.mzMaxDataStructure,
                &tandemScanIons
        );

        QVector<QVector<TandemScanIon>> tranchedTandemScanIonsVec;
        e = ParallelUtils::tranchVectorForParallelizationInOrder(
                tandemScanIons,
                pythiaParameters.chunkSize,
                0,
                &tranchedTandemScanIonsVec
        ); ree;

        for (int i = 0; i < tranchedTandemScanIonsVec.size(); i++) {
            const QVector<TandemScanIon> &tandemScansTranche = tranchedTandemScanIonsVec.at(i);
            tranchedTandemScanIons->insert(i, tandemScansTranche);
        }

        ERR_RETURN
    }

    QPair<Err, FragmentLibraryRTree*> buildRtreesParallelLogic(const FragLibIonTranche &flit) {

        ERR_INIT

        auto *rTree = new FragmentLibraryRTree();
        rTree->setKey(flit.key);
        e = rTree->init(
                flit.fragLibIonsTranche,
                flit.minMaxCharge,
                flit.ppmTolerance,
                flit.precursorExtractionWindowThomsons
                );
        // NOTE: error is handled at calling function.
        return {e, rTree};
    }

    struct RowToWrite {
        ScanNumber scanNumber = -1;
        PeptideId peptideId = -1;
        Occurrence occurrence = -1;
        int intensity = 0;
        double meanMzPM = -1.0;
        QString sequence;
        QString sequenceWithMods;
        QChar previousResidue;
        QChar postResidue;
        double mass = -1.0;
        bool isDecoy = false;
    };

    Err matchPeptideIdstoPeptides(
            const QString &pepLibFileURI,
            const QMap<ScanNumber, QVector<TallyPeptideId>> &tallyItemsByScanNumber,
            QVector<RowToWrite> *rowsToWrite
            ) {

        ERR_INIT

        PeptidesLibraryTron peptidesLibraryTron;
        e = peptidesLibraryTron.readPeptidesLib(pepLibFileURI); ree;
        qDebug() << "SDFSD" << tallyItemsByScanNumber.size();
        for (auto it = tallyItemsByScanNumber.begin(); it != tallyItemsByScanNumber.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const QVector<TallyPeptideId> &tallyPepIds = it.value();
            qDebug() << "SDFDS" << tallyPepIds.size();
            for (const TallyPeptideId &tpi : tallyPepIds) {

                e = ErrorUtils::isEqual(scanNumber, tpi.scanNumber); ree;

                Peptide peptide;
                e = peptidesLibraryTron.getPeptideById(
                        tpi.peptideId,
                        &peptide
                        ); ree;

                RowToWrite rowToWrite;
                rowToWrite.scanNumber = scanNumber;
                rowToWrite.peptideId = tpi.peptideId;
                rowToWrite.occurrence = tpi.occurrence;
                rowToWrite.intensity = tpi.intensity;
                rowToWrite.meanMzPM = tpi.meanMzPPM;
                rowToWrite.sequence = peptide.sequence;
                rowToWrite.sequenceWithMods = peptide.peptideStringWithMods();
                rowToWrite.mass = peptide.mass;
                rowToWrite.isDecoy = peptide.isDecoy;
                rowToWrite.previousResidue = peptide.previousResidue;
                rowToWrite.postResidue = peptide.postResidue;

                rowsToWrite->push_back(rowToWrite);
            }

        }

        ERR_RETURN
    }

    void deleteRTreePointers(const QMap<int, FragmentLibraryRTree*> &rTreesByKey) {
        for (FragmentLibraryRTree *rt : rTreesByKey) {
            delete rt;
        }
    }

}//namespace
Err MsFraggerTronWorkFlow::processFile(const QString &mzmLFileURI) {

    ERR_INIT

    QMap<int, QVector<TandemScanIon>> tranchedTandemScanIons;
    e = trancheTandemScanIons(
            m_pythiaParameters,
            mzmLFileURI,
            &tranchedTandemScanIons
            ); ree;

    QMap<int, FragmentLibraryRTree*> rTreesByKey;
    e = buildRTrees(
            tranchedTandemScanIons,
            &rTreesByKey
            ); ree;

    QHash<int, QVector<PeptideIdIonFraggerResult>> peptideIdIonFraggerResults;
    e = fragScanIons(
            tranchedTandemScanIons,
            rTreesByKey,
            &peptideIdIonFraggerResults
            ); ree;

    QMap<ScanNumber, QVector<TallyPeptideId>> tallyItemsByScanNumber;
    e = tallyPeptideIdsPerScan(
            peptideIdIonFraggerResults,
            &tallyItemsByScanNumber
            ); ree;

    QVector<RowToWrite> rowsToWrite;
    e = matchPeptideIdstoPeptides(
            m_pepLibUri,
            tallyItemsByScanNumber,
            &rowsToWrite
            ); ree;

    deleteRTreePointers(rTreesByKey);

    ERR_RETURN
}

namespace {

    QMap<int, QPair<MzMin, MzMax>> getEachTranchedScanIonsMzLimits(
            const QMap<int, QVector<TandemScanIon>> &tranchedTandemScanIons
    ) {

        const auto sortLogic
                = [](const TandemScanIon &l, const TandemScanIon &r) {return l.mz < r.mz;};

        QMap<int, QPair<MzMin, MzMax>> minMaxLimitsPerTranche;
        for (auto it = tranchedTandemScanIons.begin(); it != tranchedTandemScanIons.end(); it++){

            const int key = it.key();
            const QVector<TandemScanIon> &vec = it.value();

            const auto minMaxMz = std::minmax_element(vec.begin(), vec.end(), sortLogic);
            minMaxLimitsPerTranche.insert(key, {minMaxMz.first->mz, minMaxMz.second->mz});
        }

        return minMaxLimitsPerTranche;
    }

    Err buildFragLibIonTranches(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QMap<int, QPair<MzMin, MzMax>> &mzMinMaxLimitsOfTranches,
            QMap<int, FragLibIonTranche> *fragLibIonTranches
    ) {

        ERR_INIT

        FragLibraryTron fragLibraryTron(pythiaParameters);
        e = fragLibraryTron.readFragLibIons(fragLibUri); ree;
        e = fragLibraryTron.getFragLibIonTranches(
                mzMinMaxLimitsOfTranches,
                fragLibIonTranches
        ); ree;

        ERR_RETURN
    }

}//namespace
Err MsFraggerTronWorkFlow::buildRTrees(
        const QMap<int, QVector<TandemScanIon>> &tranchedTandemScanIons,
        QMap<int, FragmentLibraryRTree *> *rTreesByKey
        ) {

    ERR_INIT

    const QMap<int, QPair<MzMin, MzMax>> mzMinMaxLimitsOfTranches
            = getEachTranchedScanIonsMzLimits(tranchedTandemScanIons);

    QMap<int, FragLibIonTranche> fragLibIonTranches;
    e = buildFragLibIonTranches(
            m_pythiaParameters,
            m_fragLibUri,
            mzMinMaxLimitsOfTranches,
            &fragLibIonTranches
    ); ree;

#define RUN_PARALLEL_RTREE_LOADING
#ifdef RUN_PARALLEL_RTREE_LOADING
    QFuture<QPair<Err, FragmentLibraryRTree*>> futures = QtConcurrent::mapped(
            fragLibIonTranches,
            buildRtreesParallelLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, FragmentLibraryRTree*> &res : futures) {
        e = res.first; ree;
        rTreesByKey->insert(res.second->getKey(), res.second);
    }
#else
    for (auto it = fragLibIonTranches.begin(); it != fragLibIonTranches.end(); it++) {

        const int key = it.key();
        const FragLibIonTranche &flit = it.value();

        qDebug() << "Loading rtree key" << key;

        QPair<Err, FragmentLibraryRTree*> rTreeResult = buildRtreesParallelLogic(flit);
        e = rTreeResult.first; ree;

        rTreesByKey->insert(rTreeResult.second->getKey(), rTreeResult.second);
    }
#endif

    e = ErrorUtils::isEqual(rTreesByKey->size(), fragLibIonTranches.size()); ree;

    ERR_RETURN
}

namespace {

    struct FraggerParallelInput {
        FragmentLibraryRTree *rTree;
        QVector<TandemScanIon> tandemScanIons;
    };

    Err buildFraggerParallelInputs(
            const QMap<FraggerKey, QVector<TandemScanIon>> &tranchedTandemScanIons,
            const QMap<FraggerKey, FragmentLibraryRTree *> &rTreesByKey,
            QVector<FraggerParallelInput> *fraggerParallelInputs
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(tranchedTandemScanIons); ree;
        e = ErrorUtils::isEqual(tranchedTandemScanIons.size(), rTreesByKey.size()); ree;

        for (FraggerKey key : tranchedTandemScanIons.keys()) {

            e = ErrorUtils::isTrue(rTreesByKey.contains(key)); ree;

            FraggerParallelInput fraggerParallelInput;
            fraggerParallelInput.rTree = rTreesByKey.value(key);
            fraggerParallelInput.tandemScanIons = tranchedTandemScanIons.value(key),

            fraggerParallelInputs->push_back(fraggerParallelInput);
        }

        ERR_RETURN
    }

    QVector<PeptideIdIonFraggerResult> fraggerParallelLogic(const FraggerParallelInput &fraggerParallelInput) {

        QVector<PeptideIdIonFraggerResult> peptideIdIonResults;

        for (const TandemScanIon &tsi : fraggerParallelInput.tandemScanIons) {

            const QPair<double, double> targetWindow = {tsi.precursorTargetLowerWindow, tsi.precursorTargetUpperWindow};

            const QHash<PeptideId, MZION> peptidesTableIdsSearchResult = fraggerParallelInput.rTree->getPeptidesTableIds(
                    tsi.mz,
                    tsi.precursorTargetMz,
                    targetWindow);

            for (auto it = peptidesTableIdsSearchResult.begin(); it != peptidesTableIdsSearchResult.end(); it++) {

                const PeptideId pepId = it.key();
                const double rTreeMz = it.value();
                const double ppmMzSearched = (std::abs(rTreeMz - tsi.mz) / rTreeMz) * 1e6;

                PeptideIdIonFraggerResult res;
                res.scanNumber = tsi.scanNumber;
                res.peptideId = pepId;
                res.searchedFragIonMz = tsi.mz;
                res.intensity = tsi.intensity;
                res.ppmMzSearched = ppmMzSearched;
                peptideIdIonResults.push_back(res);
            }

        }

        return peptideIdIonResults;
    }


}//namespace
Err MsFraggerTronWorkFlow::fragScanIons(
        const QMap<FraggerKey, QVector<TandemScanIon>> &tranchedTandemScanIons,
        const QMap<FraggerKey, FragmentLibraryRTree *> &rTreesByKey,
        QHash<ScanNumber , QVector<PeptideIdIonFraggerResult>> *peptideIdIonFraggerResults
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(tranchedTandemScanIons); ree;
    e = ErrorUtils::isNotEmpty(rTreesByKey); ree;

    QVector<FraggerParallelInput> fraggerParallelInputs;
    e = buildFraggerParallelInputs(
            tranchedTandemScanIons,
            rTreesByKey,
            &fraggerParallelInputs
            ); ree;

#define RUN_PARALLEL_FRAGGER
#ifdef RUN_PARALLEL_FRAGGER
    QFuture<QVector<PeptideIdIonFraggerResult>> futures = QtConcurrent::mapped(
            fraggerParallelInputs,
            fraggerParallelLogic
            );
    futures.waitForFinished();

    for (const QVector<PeptideIdIonFraggerResult> &future : futures) {
        for (const PeptideIdIonFraggerResult &res : future) {
            (*peptideIdIonFraggerResults)[res.scanNumber].push_back(res);
        }
    }
#else
    for (const FraggerParallelInput &fpi : fraggerParallelInputs) {
        qDebug() << "RTree id:" << fpi.rTree->getKey() << "Rtree size:" << fpi.rTree->size();
        const QVector<PeptideIdIonFraggerResult> future = fraggerParallelLogic(fpi);
        for (const PeptideIdIonFraggerResult &res : future) {
            (*peptideIdIonFraggerResults)[res.scanNumber].push_back(res);
        }
    }
#endif

    ERR_RETURN
}

namespace {

    void fiterTallyPeptideIdByOccurrences(QVector<TallyPeptideId> *tallyPeptideIdsVec) {

        const int minOccurrenceCount = 3;
        const auto terminatorLogic = [minOccurrenceCount](const TallyPeptideId &tpi){
            return tpi.occurrence < minOccurrenceCount;
        };

        const auto terminator = std::remove_if(
                tallyPeptideIdsVec->begin(),
                tallyPeptideIdsVec->end(),
                terminatorLogic
                );

        tallyPeptideIdsVec->erase(terminator, tallyPeptideIdsVec->end());
    }

    void returnTopNTallyPeptideIds(QVector<TallyPeptideId> *tallyPeptideIdsVec) {

        const int occurrencesTopN = 50; //TODO consider making this settable.

        const auto occurrencesTopNLogic = [](const TallyPeptideId &l, const TallyPeptideId &r){

            if (l.occurrence == r.occurrence) {

                if (static_cast<int>(l.intensity) == static_cast<int>(r.intensity)) {
                    return l.meanMzPPM < r.meanMzPPM;
                }

                return l.intensity > r.intensity;
            }

            return l.occurrence > r.occurrence;
        };

        const int maxSize = std::min(occurrencesTopN, tallyPeptideIdsVec->size());

        std::sort(
                tallyPeptideIdsVec->begin(),
                tallyPeptideIdsVec->end(),
                occurrencesTopNLogic
                );

        tallyPeptideIdsVec->resize(maxSize);
    }

    QVector<TallyPeptideId> tallyFraggerResults(const QVector<PeptideIdIonFraggerResult> &peptideIdIonFraggerResults) {

        QHash<PeptideId, TallyPeptideId> tallyPeptideIds;

        for (const PeptideIdIonFraggerResult &pifr : peptideIdIonFraggerResults) {

            TallyPeptideId &tallyPeptideId = tallyPeptideIds[pifr.peptideId];

            if (tallyPeptideId.peptideId == -1) {
                tallyPeptideId.peptideId = pifr.peptideId;
                tallyPeptideId.scanNumber = pifr.scanNumber;
                tallyPeptideId.meanMzPPM = 0;
            }

            tallyPeptideId.meanMzPPM
                = ((tallyPeptideId.meanMzPPM * tallyPeptideId.occurrence) + pifr.ppmMzSearched) / (tallyPeptideId.occurrence + 1);

            tallyPeptideId.occurrence++;
            tallyPeptideId.intensity += pifr.intensity;
        }

        QVector<TallyPeptideId> tallyPeptideIdsVec = tallyPeptideIds.values().toVector();
        fiterTallyPeptideIdByOccurrences(&tallyPeptideIdsVec);
        returnTopNTallyPeptideIds(&tallyPeptideIdsVec);

        return tallyPeptideIdsVec;
    }


}//namespace
Err MsFraggerTronWorkFlow::tallyPeptideIdsPerScan(
        const QHash<ScanNumber, QVector<PeptideIdIonFraggerResult>> &peptideIdIonFraggerResults,
        QMap<ScanNumber, QVector<TallyPeptideId>> *tallyItemsByScanNumber
        ) {

    ERR_INIT

#define PARALLEL_PROCESS_TALLY_ITEMS
#ifdef PARALLEL_PROCESS_TALLY_ITEMS
    QFuture<QVector<TallyPeptideId>> futures = QtConcurrent::mapped(
            peptideIdIonFraggerResults,
            tallyFraggerResults
            );
    futures.waitForFinished();

    for (const QVector<TallyPeptideId> &tpi : futures) {

        if (tpi.isEmpty()) {
            continue;
        }

        (*tallyItemsByScanNumber)[tpi.front().scanNumber] = tpi;
    }
#else
    for (const QVector<PeptideIdIonFraggerResult> &pifr : peptideIdIonFraggerResults) {

        QVector<TallyPeptideId> tallyPeptideIds = tallyFraggerResults(pifr);

        if (tallyPeptideIds.isEmpty()) {
            continue;
        }

        (*tallyItemsByScanNumber)[tallyPeptideIds.front().scanNumber] = tallyPeptideIds;
    }
#endif

    ERR_RETURN
}
