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
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide
            ) {

        ERR_INIT

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msDataFilePath); ree;

        const QVector<MsScanInfo> msScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();
        const RTree rtree = loadScanInfoToRTree(msScanInfos);

        FragLibReader fragLibReader;
        e = fragLibReader.init(fragLibFilePath); ree;

        QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptide;


        e = fragLibReader.getMS2IonsTopN(
                topNMs2Ions,
                &peptideSequenceChargeKeyVsCandidatePeptide
                ); ree;

        for (auto it = peptideSequenceChargeKeyVsCandidatePeptide.begin(); it != peptideSequenceChargeKeyVsCandidatePeptide.end(); it++) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            const CandidatePeptide &candidatePeptide = it.value();

            const double mz = BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge);

            const rTreeBox queryBox(
                    rTreeCoor(mz, 0.0),
                    rTreeCoor(mz, 0.0)
            );

            std::vector<rTreePoint> rTreeSearchResult;
            rtree.query(bgi::intersects(queryBox), std::back_inserter(rTreeSearchResult));

            for (const rTreePoint &rtp : rTreeSearchResult) {
                const UniqueMsInfoScanKey &uniqueMsInfoScanKey = rtp.second;

                (*uniqueInfoScanKeyVsCandidatePeptide)[uniqueMsInfoScanKey].insert(candidatePeptide.peptideStringWithMods, candidatePeptide);
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
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptide = nullptr;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        qDebug() << "Processing" << ppi.uniqueMsInfoScanKey;

        MsFrameScoretron msFrameScoretron;
        e = msFrameScoretron.init(
                ppi.pythiaParameters,
                *ppi.featureFinderHills,
                *ppi.peptideStringWithModsVsCandidatePeptide,
                ppi.scanNumberVsScanTime,
                ppi.iRTReCalFilePath
                ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        qDebug() << "Processed" << ppi.uniqueMsInfoScanKey << "Count" << scoredCandidates.size() << et.elapsed() << "mSec";

        return {e, scoredCandidates};
    }


}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    QMap<UniqueMsInfoScanKey, QVector<FeatureFinderHill>> uniqueInfoScanKeyVsFeatureFinderHills;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    e = buildUniqueInfoScanKeyVsFeatureFinderHills(
            msDataFilePath,
            &uniqueInfoScanKeyVsFeatureFinderHills,
            &scanNumberVsScanTime
            ); ree;

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptide;
    e = buildLibrary(
            msDataFilePath,
            m_fragLibUri,
            m_pythiaParameters.topNMs2Ions,
            &uniqueInfoScanKeyVsCandidatePeptide
            ); ree;

    QVector<ParallelProcessingInput> parallelProcessingInputs;
    for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueInfoScanKeyVsCandidatePeptide.keys()) {
        ParallelProcessingInput ppi;
        ppi.uniqueMsInfoScanKey = uniqueMsInfoScanKey;
        ppi.scanNumberVsScanTime = scanNumberVsScanTime;
        ppi.iRTReCalFilePath = m_iRTReCalFilePath;
        ppi.pythiaParameters = m_pythiaParameters;
        ppi.featureFinderHills = &uniqueInfoScanKeyVsFeatureFinderHills[uniqueMsInfoScanKey];
        ppi.peptideStringWithModsVsCandidatePeptide = &uniqueInfoScanKeyVsCandidatePeptide[uniqueMsInfoScanKey];

        parallelProcessingInputs.push_back(ppi);
    }

    QVector<ScoredCandidate> combinedResults;

//#define PROCESS_PARALLEL
#ifdef PROCESS_PARALLEL
    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProciessingLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<ScoredCandidate>> &res : futures) {

        e = res.first; ree;
        combinedResults.append(res.second);
    }
#else
    for (const ParallelProcessingInput &inp : parallelProcessingInputs) {
        const QPair<Err, QVector<ScoredCandidate>> res = parallelProciessingLogic(inp);
        e = res.first; ree;
        combinedResults.append(res.second);
    }
#endif

    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
    e = ParquetReader::write(combinedResults, resultsFilePath); ree;

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

    Err loadHillsToScanNumberVsScanPoints(
            const QVector<FeatureFinderHill> &featureFinderHills,
            QMap<ScanNumber, ScanPoints> *scanNumberScanPointsHills
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

        for (const FeatureFinderHill &ffh : featureFinderHills) {

            const double mzMean = ffh.mzMean();
            const QVector<ScanNumber> &scanNumbers = ffh.scanNumbers();
            const QVector<double> &intensities = ffh.intensities();

            e = ErrorUtils::isEqual(scanNumbers.size(), intensities.size()); ree

            for (int i = 0; i < scanNumbers.size(); i++) {

                const ScanNumber scanNumber = scanNumbers.at(i);
                const double intensity = intensities.at(i);

                if (MathUtils::tZero(intensity)) {
                    continue;
                }

                (*scanNumberScanPointsHills)[scanNumber].push_back({mzMean, intensity});
            }
        }

        ERR_RETURN
    }

    Err buildSmoothedScanNumberVsScanPoints(
            const QVector<FeatureFinderHill> &featureFinderHills,
            const FeatureFinderParameters &ffParams,
            double mzMax,
            double ppmTol,
            QMap<ScanNumber, ScanPoints> *smoothedScanNumberVsScanPoints
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

        QMap<ScanNumber, ScanPoints> scanNumberScanPointsHills;
        e = loadHillsToScanNumberVsScanPoints(featureFinderHills, &scanNumberScanPointsHills); ree;

        MsFrame msFrame;
        msFrame.init(
                "KalliopeAndBellatrix",
                scanNumberScanPointsHills,
                {-1.0, -1.0},
                {{-1, -1.0}}
                ); ree;

        e = msFrame.deisotopeMsFrame(ppmTol); ree;
        e = msFrame.smoothFrame(
                ffParams.filterLength,
                ffParams.sigma,
                ffParams.smoothCount,
                mzMax
                ); ree;

        *smoothedScanNumberVsScanPoints = msFrame.scanNumberVsScanPoints();

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

        QVector<FeatureFinderHill> featureFinderHills;
        e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills); rree;

        QMap<ScanNumber, ScanPoints> smoothedScanNumberVsScanPoints;
        e = buildSmoothedScanNumberVsScanPoints(
                featureFinderHills,
                ffParams,
                pythiaParameters.mzMaxDataStructure,
                pythiaParameters.ms2ExtractionWidthPPM,
                &smoothedScanNumberVsScanPoints
                ); rree;

        e = featureFinderHillBuilder.buildHills(smoothedScanNumberVsScanPoints); rree;
        e = featureFinderHillBuilder.refineHills(false); rree;

        QVector<FeatureFinderHill> smoothedFeatureFinderHills;
        e = featureFinderHillBuilder.featureFinderHills(&smoothedFeatureFinderHills); rree;

        return {e, {frame.first, smoothedFeatureFinderHills}};
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



