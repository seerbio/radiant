//
// Created by anichols on 2/19/23.
//

#include "MsFraggerTronWorkFlow.h"

#include "DeisotoperTandem.h"
#include "MsScansDenoiseTron.h"
#include "FragLibraryTron.h"
#include "MsFraggerTronResultsReader.h"
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

    Err trancheTandemScanIons(
            const PythiaParameters &pythiaParameters,
            const QVector<TandemScanIon> &tandemScanIons,
            QMap<Tranche, QVector<TandemScanIon>> *tranchedTandemScanIons
            ) {

        ERR_INIT

        const int trancheBuffer = 0;

        QVector<QVector<TandemScanIon>> tranchedTandemScanIonsVec;
        e = ParallelUtils::tranchVectorForParallelizationInOrder(
                tandemScanIons,
                pythiaParameters.trancheSize,
                trancheBuffer,
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

    void deleteRTreePointers(const QMap<int, FragmentLibraryRTree*> &rTreesByKey) {
        for (FragmentLibraryRTree *rt : rTreesByKey) {
            delete rt;
        }
    }

}//namespace
Err MsFraggerTronWorkFlow::processFile(
        const QString &mzmLFileURI,
        QString *psmResultsFilePath
        ) {

    ERR_INIT

    MsReaderMzML mzMLReader;
    e = mzMLReader.openFile(mzmLFileURI); ree;

    QVector<TandemScanIon> tandemScanIons;
    e = mzMLReader.tandemScanIons(&tandemScanIons); ree;

    e = preProcessScans(&tandemScanIons); ree;

    QMap<ScanNumber, QVector<TallyPeptideId>> tallyItemsByScanNumber;
    e = processScanIons(
            tandemScanIons,
            &tallyItemsByScanNumber
            ); ree;

    *psmResultsFilePath = mzmLFileURI + S_GLOBAL_SETTINGS.DOT_PSM + S_GLOBAL_SETTINGS.DOT_CSV;

    e = writePSMsToFile(
            *psmResultsFilePath,
            tallyItemsByScanNumber
            ); ree;

    ERR_RETURN
}

Err MsFraggerTronWorkFlow::processScanIons(
        const QVector<TandemScanIon> &tandemScanIons,
        QMap<ScanNumber, QVector<TallyPeptideId>> *tallyItemsByScanNumber
        ) {

    ERR_INIT

    QMap<Tranche, QVector<TandemScanIon>> tranchedTandemScanIons;
    e = trancheTandemScanIons(
            m_pythiaParameters,
            tandemScanIons,
            &tranchedTandemScanIons
    ); ree;

    QMap<int, FragmentLibraryRTree*> rTreesByKey;
    e = buildRTrees(
            tranchedTandemScanIons,
            &rTreesByKey
    ); ree;

    QHash<ScanNumber, QVector<PeptideIdIonFraggerResult>> peptideIdIonFraggerResults;
    e = fragScanIons(
            tranchedTandemScanIons,
            rTreesByKey,
            &peptideIdIonFraggerResults
    ); ree;

    deleteRTreePointers(rTreesByKey);

    e = tallyPeptideIdsPerScan(
            peptideIdIonFraggerResults,
            tallyItemsByScanNumber
    ); ree;

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

    struct DeisotopeParallelInput {
        QMap<ScanNumber, ScanPoints> scanPointsFrame;
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        double ppmTol = -1.0;
    };

    QVector<DeisotopeParallelInput> buildDisotopeParallelInputs(
            const QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &diaFrames,
            double ppmTol
            ) {

        QVector<DeisotopeParallelInput> deisotopeParallelInputs;
        for(auto it = diaFrames.begin(); it != diaFrames.end(); it++) {

            DeisotopeParallelInput dipi;
            dipi.scanPointsFrame = it.value();
            dipi.uniqueMsInfoScanKey = it.key();
            dipi.ppmTol = ppmTol;

            deisotopeParallelInputs.push_back(dipi);
        }

        return deisotopeParallelInputs;
    }

    QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> deisotopeTargetMzScansFrame (
            const DeisotopeParallelInput &deisotopeParallelInput
            ) {

        const QMap<ScanNumber, ScanPoints> &sps = deisotopeParallelInput.scanPointsFrame;

        QMap<ScanNumber, ScanPoints> deisotopedScanFrame;
        for (auto it = sps.begin(); it != sps.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const ScanPoints &sp = it.value();

            const ScanPoints deisotopedScanPoints
                = DeisotoperTandem::deisotopeTandemScan(sp, deisotopeParallelInput.ppmTol);

            deisotopedScanFrame.insert(scanNumber, deisotopedScanPoints);
        }

        return {deisotopeParallelInput.uniqueMsInfoScanKey, deisotopedScanFrame};
    }

    QList<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> deisotopeScanFrames(
            const QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &diaFrames,
            const PythiaParameters &pythiaParameters
            ) {

        const QVector<DeisotopeParallelInput> deisotopeParallelInputs = buildDisotopeParallelInputs(
                diaFrames,
                pythiaParameters.ms2ExtractionWidthPPM
        );

        QFuture<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> futures = QtConcurrent::mapped(
                deisotopeParallelInputs,
                deisotopeTargetMzScansFrame
        );
        futures.waitForFinished();

        return futures.results();
    }

    Err denoiseScanFramesParallelLogic(
            const QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &scanFrame,
            const PythiaParameters &pythiaParameters,
            QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *denoisedScanFrame
            ) {

        ERR_INIT

        MsScansDenoiseTron msScansDenoiseTron;

        e = msScansDenoiseTron.init(pythiaParameters); ree;

        QMap<ScanNumber, ScanPoints> scansFrameDenoised;
        e = msScansDenoiseTron.denoiseScansFrame(
                scanFrame.second,
                &scansFrameDenoised
                ); ree;

        *denoisedScanFrame = {scanFrame.first, scansFrameDenoised};

        ERR_RETURN
    }

    Err denoiseScanFrames(
            const QList<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> &scanFrames,
            const PythiaParameters &pythiaParameters,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *denoisedScanFrames
            ) {

        ERR_INIT

        for (const QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &scanFrame : scanFrames) {

            QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> targetDenoisedScanFrame;
            e = denoiseScanFramesParallelLogic(
                    scanFrame,
                    pythiaParameters,
                    &targetDenoisedScanFrame
                    ); ree;

            denoisedScanFrames->insert(
                    targetDenoisedScanFrame.first,
                    targetDenoisedScanFrame.second
                    );
        }

        ERR_RETURN
    }

    QHash<UniqueMsInfoScanKey, TandemScanIon> buildTandemScanReconstructionHash(
            const QVector<TandemScanIon> &tandemScanIons
    ) {

        QHash<UniqueMsInfoScanKey, TandemScanIon> reconMap;

        for (const TandemScanIon &tsi : tandemScanIons) {

            const UniqueMsInfoScanKey uniqueMsInfoScanKey = tsi.targetScanKey();

            if (reconMap.value(uniqueMsInfoScanKey).scanNumber != -1) {
                continue;
            }

            reconMap.insert(uniqueMsInfoScanKey, tsi);
        }

        return reconMap;
    }

    Err rejoinTandemScanInfoWithProcessedScanPoints(
            const QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &denoisedScanFrames,
            QVector<TandemScanIon> *tandemScanIons
            ) {

        ERR_INIT

        const QHash<UniqueMsInfoScanKey, TandemScanIon> tandemScanIonsReconstructionHash
                = buildTandemScanReconstructionHash(*tandemScanIons);

        tandemScanIons->clear();

        for (auto it = denoisedScanFrames.begin(); it != denoisedScanFrames.end(); it++) {

            const UniqueMsInfoScanKey &uniqueMsInfoScanKey = it.key();
            const QMap<ScanNumber, ScanPoints> &dnsf = it.value();

            const TandemScanIon &reconHashRef = tandemScanIonsReconstructionHash.value(uniqueMsInfoScanKey);
            e = ErrorUtils::isNotEqual(reconHashRef.scanNumber, -1); ree;

            for (auto itt = dnsf.begin(); itt != dnsf.end(); itt++) {

                const ScanNumber scanNumber = itt.key();
                const ScanPoints &scanPoints = itt.value();

                for (const ScanPoint &sp: scanPoints) {

                    TandemScanIon newTandemScanIon;
                    newTandemScanIon.scanNumber = scanNumber;
                    newTandemScanIon.mz = sp.x();
                    newTandemScanIon.intensity = sp.y();
                    newTandemScanIon.precursorTargetMz = reconHashRef.precursorTargetMz ;
                    newTandemScanIon.precursorTargetLowerWindow = reconHashRef.precursorTargetLowerWindow;
                    newTandemScanIon.precursorTargetUpperWindow = reconHashRef.precursorTargetUpperWindow;

                    tandemScanIons->push_back(newTandemScanIon);
                }
            }
        }

        ERR_RETURN
    }

    void sortTandemScanIons(QVector<TandemScanIon> *tandemScanIons) {

        QMap<NominalMzMass, QVector<TandemScanIon>> tandemScanIonsMapped;
        for (const TandemScanIon &tsi : *tandemScanIons) {
            tandemScanIonsMapped[static_cast<int>(std::round(tsi.mz))].push_back(tsi);
        }

        MsReaderBase::sortTandemScanIonsInChunks(&tandemScanIonsMapped);

        tandemScanIons->clear();

        for (const QVector<TandemScanIon> &vec : tandemScanIonsMapped) {
            tandemScanIons->append(vec);
        }
    }

}//namespace
Err MsFraggerTronWorkFlow::preProcessScans(QVector<TandemScanIon> *tandemScanIons) {

    ERR_INIT
    qDebug() << "Before Scans Preprocessing size" << tandemScanIons->size();

    filterTandemScanIonsByMz(
            m_pythiaParameters.mzMinDataStructure,
            m_pythiaParameters.mzMaxDataStructure,
            tandemScanIons
    );

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaFrames;
    e = MsReaderMzML::sortDIATandemScansByMzTarget(
            *tandemScanIons,
            &diaFrames
    ); ree;

    const QList<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> deisotopedScanFrames = deisotopeScanFrames(
            diaFrames,
            m_pythiaParameters
            );

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> denoisedScanFrames;
    e = denoiseScanFrames(
            deisotopedScanFrames,
            m_pythiaParameters,
            &denoisedScanFrames
            ); ree;

    e = rejoinTandemScanInfoWithProcessedScanPoints(
            denoisedScanFrames,
            tandemScanIons
            ); ree;

    sortTandemScanIons(tandemScanIons);

    qDebug() << "After Scans Preprocessing size" << tandemScanIons->size();
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
        const QMap<Tranche, QVector<TandemScanIon>> &tranchedTandemScanIons,
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
                res.foundTheoFragIonMz = rTreeMz;
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

    void returnTopNTallyPeptideIds(
            int occurrencesTopN,
            QVector<TallyPeptideId> *tallyPeptideIdsVec
            ) {

        const auto occurrencesTopNLogic = [](const TallyPeptideId &l, const TallyPeptideId &r){

            if (l.occurrence == r.occurrence) {

                if (static_cast<int>(l.intensityTotal) == static_cast<int>(r.intensityTotal)) {
                    return l.meanMzPPM < r.meanMzPPM;
                }

                return l.intensityTotal > r.intensityTotal;
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

    QVector<TallyPeptideId> tallyFraggerResults(
            const QVector<PeptideIdIonFraggerResult> &peptideIdIonFraggerResults,
            int occurrencesTopN,
            bool saveMzVals
            ) {

        QHash<PeptideId, TallyPeptideId> tallyPeptideIds;

        for (const PeptideIdIonFraggerResult &pifr : peptideIdIonFraggerResults) {

            TallyPeptideId &tallyPeptideId = tallyPeptideIds[pifr.peptideId];

            if (tallyPeptideId.peptideId == -1) {
                tallyPeptideId.peptideId = pifr.peptideId;
                tallyPeptideId.scanNumber = pifr.scanNumber;
                tallyPeptideId.meanMzPPM = 0;
            }

            tallyPeptideId.meanMzPPM
                = ((tallyPeptideId.meanMzPPM * tallyPeptideId.occurrence) + pifr.ppmMzSearched)
                        / (tallyPeptideId.occurrence + 1);

            tallyPeptideId.occurrence++;
            tallyPeptideId.intensityTotal += static_cast<int>(pifr.intensity);
            if (saveMzVals) {
                tallyPeptideId.scanIonMZs.push_back(pifr.searchedFragIonMz);
                tallyPeptideId.theoFragMZs.push_back(pifr.foundTheoFragIonMz);
            }
        }

        QVector<TallyPeptideId> tallyPeptideIdsVec = tallyPeptideIds.values().toVector();
        fiterTallyPeptideIdByOccurrences(&tallyPeptideIdsVec);
        returnTopNTallyPeptideIds(occurrencesTopN, &tallyPeptideIdsVec);

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

    const auto tallyFraggerSaveBinder = std::bind(
            tallyFraggerResults,
            std::placeholders::_1,
            m_pythiaParameters.returnPSMTopN,
            true
            );

    QFuture<QVector<TallyPeptideId>> futures = QtConcurrent::mapped(
            peptideIdIonFraggerResults,
            tallyFraggerSaveBinder
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

namespace {

    Err matchPeptideIdstoPeptides(
            const QString &pepLibFileURI,
            const QMap<ScanNumber, QVector<TallyPeptideId>> &tallyItemsByScanNumber,
            QVector<RowToWrite> *rowsToWrite
    ) {

        ERR_INIT

        PeptidesLibraryTron peptidesLibraryTron;
        e = peptidesLibraryTron.readPeptidesLib(pepLibFileURI); ree;

        for (auto it = tallyItemsByScanNumber.begin(); it != tallyItemsByScanNumber.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const QVector<TallyPeptideId> &tallyPepIds = it.value();

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
                rowToWrite.intensityTotal = tpi.intensityTotal;
                rowToWrite.meanMzPPM = tpi.meanMzPPM;
                rowToWrite.sequence = peptide.sequence;
                rowToWrite.sequenceWithMods = peptide.peptideStringWithMods();
                rowToWrite.mass = peptide.mass;
                rowToWrite.isDecoy = peptide.isDecoy;
                rowToWrite.previousResidue = peptide.previousResidue;
                rowToWrite.postResidue = peptide.postResidue;
                rowToWrite.scanIonMZs = tpi.scanIonMZs;
                rowToWrite.theoFragIonMZs = tpi.theoFragMZs;
                rowsToWrite->push_back(rowToWrite);
            }

        }

        ERR_RETURN
    }

}//namespace
Err MsFraggerTronWorkFlow::writePSMsToFile(
        const QString &outputFilePath,
        const QMap<ScanNumber, QVector<TallyPeptideId>> &tallyItemsByScanNumber
        ) {

    ERR_INIT

    QVector<RowToWrite> rowsToWrite;
    e = matchPeptideIdstoPeptides(
            m_pepLibUri,
            tallyItemsByScanNumber,
            &rowsToWrite
    ); ree;

    MsFraggerTronResultsReader::writeToCsv(
            outputFilePath,
            &rowsToWrite
    );

    ERR_RETURN
}
