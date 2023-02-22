//
// Created by anichols on 2/19/23.
//

#include "MsFraggerTronWorkFlow.h"

#include "FragLibraryTron.h"
#include "MsReaderMzML.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>
#include <QFuture>


Err MsFraggerTronWorkFlow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

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



    deleteRTreePointers(rTreesByKey);

    ERR_RETURN
}

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
