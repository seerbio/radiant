//
// Created by anichols on 2/19/23.
//

#include "MsFraggerTronWorkFlow.h"

#include "FragLibraryTron.h"
#include "MsReaderMzML.h"
#include "ParallelUtils.h"


MsFraggerTronWorkFlow::~MsFraggerTronWorkFlow() {
    delete m_fragLibraryTron;
}

Err MsFraggerTronWorkFlow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri) {

    ERR_INIT

    m_pythiaParameters = pythiaParameters;

    m_fragLibraryTron = new FragLibraryTron(m_pythiaParameters);
    e = m_fragLibraryTron->readFragLibIons(fragLibUri); ree;

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

    QVector<QPair<MzMin, MzMax>> getEachTranchedScanIonsMzLimits(
            const QVector<QVector<TandemScanIon>> &tranchedTandemScanIons
            ) {

        const auto sortLogic
            = [](const TandemScanIon &l, const TandemScanIon &r) {return l.mz < r.mz;};

        QVector<QPair<MzMin, MzMax>> minMaxLimitsPerTranche;
        for (const QVector<TandemScanIon> &vec : tranchedTandemScanIons){
            const auto minMaxMz = std::minmax_element(vec.begin(), vec.end(), sortLogic);
            minMaxLimitsPerTranche.push_back({minMaxMz.first->mz, minMaxMz.second->mz});
        }

        return minMaxLimitsPerTranche;
    }

}//namespace
Err MsFraggerTronWorkFlow::processFile(const QString &mzmLFileURI) {

    ERR_INIT

    QVector<QVector<TandemScanIon>> tranchedTandemScanIons;
    {
        MsReaderMzML mzMLReader;
        e = mzMLReader.openFile(mzmLFileURI); ree;

        QVector<TandemScanIon> tandemScanIons = mzMLReader.tandemScanIons();
        filterTandemScanIonsByMz(
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
                &tandemScanIons
                );

        e = ParallelUtils::tranchVectorForParallelizationInOrder(
                tandemScanIons,
                m_pythiaParameters.chunkSize,
                0,
                &tranchedTandemScanIons
                ); ree;
    }

    const QVector<QPair<MzMin, MzMax>> mzMinMaxLimitsOfTranches
            = getEachTranchedScanIonsMzLimits(tranchedTandemScanIons);

    QVector<FragLibIonTranche> fragLibIonTranches;
    e = m_fragLibraryTron->getFragLibIonTranches(
            mzMinMaxLimitsOfTranches,
            &fragLibIonTranches
            ); ree;


    ERR_RETURN
}
