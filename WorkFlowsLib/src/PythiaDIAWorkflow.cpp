//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"

#include <QtConcurrent/QtConcurrent>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
using rTreeBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeBox, QString> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri
        ) {

    ERR_INIT

    pythiaParameters.print();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    ERR_RETURN
}

Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &iRTReCalFilePath
) {

    ERR_INIT

    e = init(
            pythiaParameters,
            fragLibUri
            ); ree;

    e = ErrorUtils::fileExists(iRTReCalFilePath); ree;

    m_iRTReCalFilePath = iRTReCalFilePath;

    ERR_RETURN
}


namespace {

    RTree loadScanInfoToRTree(const QVector<MsScanInfo> &msScanInfos) {

        ERR_INIT

        std::vector<rTreePoint> cloudLoader;

        std::transform(
                msScanInfos.begin(),
                msScanInfos.end(),
                std::back_inserter(cloudLoader),
                [](const MsScanInfo &msScanInfo){

                    rTreeBox box(
                            {msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower, 0.0},
                            {msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper, 0.0}
                    );

                    return rTreePoint(box, msScanInfo.targetScanKey());
                }
            );

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    Err buildLibrary(
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            int topNMs2Ions,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *uniqueInfoScanKeyVsMs2Ions
            ) {

        ERR_INIT

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msDataFilePath); ree;

        const QVector<MsScanInfo> msScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();
        const RTree rtree = loadScanInfoToRTree(msScanInfos);

        FragLibReader fragLibReader;
        e = fragLibReader.init(fragLibFilePath); ree;

        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> peptideSequenceChargeKeyVsMS2Ions;
        QMap<PeptideSequenceChargeKey, bool> peptideSequenceChargeKeyVsIsDecoy;
        QMap<PeptideSequenceChargeKey, double> peptideSequenceChargeKeyVsMass;
        QMap<PeptideSequenceChargeKey, double> peptideSequenceChargeKeyVsIR;

        e = fragLibReader.getMS2IonsTopN(
                topNMs2Ions,
                &peptideSequenceChargeKeyVsMS2Ions,
                &peptideSequenceChargeKeyVsIsDecoy,
                &peptideSequenceChargeKeyVsMass,
                &peptideSequenceChargeKeyVsIR
                ); ree;

        for (auto it = peptideSequenceChargeKeyVsMS2Ions.begin(); it != peptideSequenceChargeKeyVsMS2Ions.end(); it++) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            const QVector<MS2Ion> &ms2IonsTopN = it.value();

            PeptideStringWithMods peptideStringWithMods;
            Charge charge;
            e = FragLibReader::peptideStringWithModsFromPeptideSequenceChargeKey(
                    peptideSequenceChargeKey,
                    &peptideStringWithMods,
                    &charge
                    ); ree;

            const double mass = peptideSequenceChargeKeyVsMass.value(peptideSequenceChargeKey);
            const double mz = BiophysicalCalcs::calculateThomsonFromMass(mass, charge);

            const rTreeBox queryBox(
                    rTreeCoor(mz, 0.0),
                    rTreeCoor(mz, 0.0)
            );

            std::vector<rTreePoint> rTreeSearchResult;
            rtree.query(bgi::intersects(queryBox), std::back_inserter(rTreeSearchResult));

            for (const rTreePoint &rtp : rTreeSearchResult) {
                const UniqueMsInfoScanKey &uniqueMsInfoScanKey = rtp.second;

                (*uniqueInfoScanKeyVsMs2Ions)[uniqueMsInfoScanKey].insert(peptideStringWithMods, ms2IonsTopN);
            }

        }

        ERR_RETURN
    }

    struct ParallelProcessingInput {
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        PythiaParameters pythiaParameters;
        QString iRTReCalFilePath;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        QVector<FeatureFinderHill> *featureFinderHills = nullptr;
        QMap<PeptideStringWithMods, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions = nullptr;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        MsFrameScoretron msFrameScoretron;
        e = msFrameScoretron.init(
                ppi.pythiaParameters,
                *ppi.featureFinderHills,
                *ppi.peptideStringWithModsVsMS2Ions,
                ppi.scanNumberVsScanTime,
                ppi.iRTReCalFilePath
                ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        return {e, scoredCandidates};
    }


}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    const QString &msDataFilePathRecalibrated = msDataFilePath; //drewholio remove this path but keep var
    qDebug() << msDataFilePathRecalibrated;

    QMap<UniqueMsInfoScanKey, QVector<FeatureFinderHill>> uniqueInfoScanKeyVsFeatureFinderHills;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    e = buildUniqueInfoScanKeyVsFeatureFinderHills(
            msDataFilePathRecalibrated,
            &uniqueInfoScanKeyVsFeatureFinderHills,
            &scanNumberVsScanTime
            ); ree;

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> uniqueInfoScanKeyVsMs2Ions;
    e = buildLibrary(
            msDataFilePathRecalibrated,
            m_fragLibUri,
            m_pythiaParameters.topNMs2Ions,
            &uniqueInfoScanKeyVsMs2Ions
            ); ree;

    QVector<ParallelProcessingInput> parallelProcessingInputs;
    for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueInfoScanKeyVsMs2Ions.keys()) {
        ParallelProcessingInput ppi;
        ppi.uniqueMsInfoScanKey = uniqueMsInfoScanKey;
        ppi.scanNumberVsScanTime = scanNumberVsScanTime;
        ppi.iRTReCalFilePath = m_iRTReCalFilePath;
        ppi.pythiaParameters = m_pythiaParameters;
        ppi.featureFinderHills = &uniqueInfoScanKeyVsFeatureFinderHills[uniqueMsInfoScanKey];
        ppi.peptideStringWithModsVsMS2Ions = &uniqueInfoScanKeyVsMs2Ions[uniqueMsInfoScanKey];

        parallelProcessingInputs.push_back(ppi);
    }

    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProciessingLogic
    );
    futures.waitForFinished();



    ERR_RETURN
}

namespace {

    Err separateTargetWindows(
            const QString &msFilePath,
            QVector<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> *frames,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
    ) {

        ERR_INIT

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msFilePath); ree;

        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrame;
        e = msReaderParquet.collateTandemPrecursorTargetsDIA(&diaTargetFrame); ree

        *frames = ParallelUtils::convertMapToVectorPairs(diaTargetFrame);

        const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderParquet.getMsScanInfos();
        for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {
            scanNumberVsScanTime->insert(it.key(), it.value().scanTime);
        }

        ERR_RETURN
    }

    QPair<Err, QPair<UniqueMsInfoScanKey, QVector<FeatureFinderHill>>> parallelFeatureFind(
            const QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &frame,
            const PythiaParameters &pythiaParameters
    ) {

        ERR_INIT

        FeatureFinderParameters ffParams;
        ffParams.tolerancePPM = pythiaParameters.ms2ExtractionWidthPPM;
        ffParams.skipScanCount = pythiaParameters.skipScanCount;
        ffParams.minScanCount = pythiaParameters.minScanCount;
        ffParams.filterLength = pythiaParameters.filterLength;
        ffParams.smoothCount = pythiaParameters.smoothCount;
        ffParams.sigma = pythiaParameters.sigma;
        ffParams.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        FeatureFinderHillBuilder featureFinderHillBuilder;
        e = featureFinderHillBuilder.init(ffParams); rree;
        e = featureFinderHillBuilder.buildHills(frame.second); rree;
//        e = featureFinderHillBuilder.refineHills(true); ree;

        QVector<FeatureFinderHill> featureFinderHills;
        e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills); rree;

        return {e, {frame.first, featureFinderHills}};
    }

}//namespace
Err PythiaDIAWorkflow::buildUniqueInfoScanKeyVsFeatureFinderHills(
        const QString &msDataFilePath,
        QMap<UniqueMsInfoScanKey, QVector<FeatureFinderHill>> *uniqueInfoScanKeyVsFeatureFinderHills,
        QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
        ) {

    ERR_INIT

    QVector<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> diaTargetFrames;
    e = separateTargetWindows(
            msDataFilePath,
            &diaTargetFrames,
            scanNumberVsScanTime
    ); ree;

#define RUN_PYTHIA_PARALLEL
#ifdef RUN_PYTHIA_PARALLEL

    const auto featureFinderLogicBinder = std::bind(
            parallelFeatureFind,
            std::placeholders::_1,
            m_pythiaParameters
    );

    QFuture<QPair<Err, QPair<UniqueMsInfoScanKey, QVector<FeatureFinderHill>>>> futures = QtConcurrent::mapped(
            diaTargetFrames,
            featureFinderLogicBinder
    );
    futures.waitForFinished();

    for (const QPair<Err, QPair<UniqueMsInfoScanKey, QVector<FeatureFinderHill>>> &res: futures) {

        e = res.first; ree;
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey = res.second.first;
        const QVector<FeatureFinderHill> &featureFinderHills = res.second.second;

        uniqueInfoScanKeyVsFeatureFinderHills->insert(uniqueMsInfoScanKey, featureFinderHills);
    }

#else
    for (const QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &targetFrame : diaTargetFrames) {
        QPair<Err, QPair<UniqueMsInfoScanKey, QVector<FeatureFinderHill>>> result = parallelFeatureFind(targetFrame, m_pythiaParameters); ree;
    }
#endif

    ERR_RETURN
}



