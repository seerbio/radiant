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
            double selectionListFraction,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide = nullptr
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

        QMap<Index, bool> selectionList;
        if (selectionListFraction > 0) {
            const int selectionCount
                = static_cast<int>(std::round(peptideSequenceChargeKeyVsCandidatePeptide.size() * selectionListFraction));

            selectionList = MathUtils::generateRandomSelectionList(peptideSequenceChargeKeyVsCandidatePeptide.size(), selectionCount);
        }

        int counter = 0;
        for (const CandidatePeptide &candidatePeptide : peptideSequenceChargeKeyVsCandidatePeptide) {

            if (!selectionList.isEmpty() && !selectionList.value(counter++)) {
                continue;
            }

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

    int calculateuniqueInfoScanKeyVsCandidatePeptideSize(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide
            ){

        const auto accumulateLogic = [](int sum, const QMap<PeptideStringWithMods, CandidatePeptide> &c){return sum + c.size();};

        return std::accumulate(
                uniqueInfoScanKeyVsCandidatePeptide.begin(),
                uniqueInfoScanKeyVsCandidatePeptide.end(),
                0,
                accumulateLogic
                );
    }

    struct ParallelProcessingInput {
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        PythiaParameters pythiaParameters;
        QString iRTReCalFilePath;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
        QMap<ScanNumber , ScanPoints> *scanPoints = nullptr;
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptide = nullptr;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        qDebug() << "Processing" << ppi.uniqueMsInfoScanKey;

        MsFrameScoretron msFrameScoretron;

        e = msFrameScoretron.init(
                ppi.uniqueMsInfoScanKey,
                ppi.pythiaParameters,
                *ppi.scanPoints,
                ppi.scanNumberVsScanTimeMS1,
                *ppi.peptideStringWithModsVsCandidatePeptide,
                ppi.scanNumberVsScanTime
        ); rree;

//        e = msFrameScoretron.init(
//                ppi.uniqueMsInfoScanKey,
//                ppi.pythiaParameters,
//                *ppi.scanPoints,
//                ppi.scanNumberVsScanTimeMS1,
//                *ppi.peptideStringWithModsVsCandidatePeptide,
//                ppi.scanNumberVsScanTime,
//                ppi.iRTReCalFilePath
//                ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        qDebug() << "Processed" << ppi.uniqueMsInfoScanKey << "Count" << scoredCandidates.size() << et.elapsed() << "mSec";

        return {e, scoredCandidates};
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    const double calibrationSelectionFraction = 0.01;
    const int minTopNMs2Ions = 6;
    const int topNMs2IonsCalibration = std::max(
            minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
            );

    qDebug() << "Using top:" << topNMs2IonsCalibration << "fragments";

    QVector<ScoredCandidate> scoredCandidatesCalibration;
    e = extractTargetDecoyData(
            msDataFilePath,
            topNMs2IonsCalibration,
            calibrationSelectionFraction,
            &scoredCandidatesCalibration
            ); ree;




    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
    e = ParquetReader::write(scoredCandidatesCalibration, resultsFilePath); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::extractTargetDecoyData(
        const QString &msDataFilePath,
        int topNMs2Ions,
        double selectionFraction,
        QVector<ScoredCandidate> *combinedResults
        ) {

    ERR_INIT

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
    e = buildLibrary(
            msDataFilePath,
            m_fragLibUri,
            topNMs2Ions,
            selectionFraction,
            &uniqueInfoScanKeyVsCandidatePeptideCalibration
    ); ree;
    qDebug() << "Calibration Size" << calculateuniqueInfoScanKeyVsCandidatePeptideSize(uniqueInfoScanKeyVsCandidatePeptideCalibration);

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> uniqueInfoScanKeyVsScanPoints;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = buildUniqueMsInfoScanKeyVsScanPoints(
            msDataFilePath,
            &uniqueInfoScanKeyVsScanPoints,
            &scanNumberVsScanTime,
            &scanNumberVsScanTimeMS1
    ); ree;

    QVector<ParallelProcessingInput> parallelProcessingInputs;
    for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueInfoScanKeyVsCandidatePeptideCalibration.keys()) {
        ParallelProcessingInput ppi;
        ppi.uniqueMsInfoScanKey = uniqueMsInfoScanKey;
        ppi.scanNumberVsScanTime = scanNumberVsScanTime;
        ppi.scanNumberVsScanTimeMS1 = scanNumberVsScanTimeMS1;
        ppi.iRTReCalFilePath = m_iRTReCalFilePath;
        ppi.pythiaParameters = m_pythiaParameters;
        ppi.scanPoints = &uniqueInfoScanKeyVsScanPoints[uniqueMsInfoScanKey];
        ppi.peptideStringWithModsVsCandidatePeptide = &uniqueInfoScanKeyVsCandidatePeptideCalibration[uniqueMsInfoScanKey];

        parallelProcessingInputs.push_back(ppi);
    }

#define PROCESS_PARALLEL
#ifdef PROCESS_PARALLEL
    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProciessingLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<ScoredCandidate>> &res : futures) {

        e = res.first; ree;
        combinedResults->append(res.second);
    }
#else
    for (const ParallelProcessingInput &inp : parallelProcessingInputs) {
        const QPair<Err, QVector<ScoredCandidate>> res = parallelProciessingLogic(inp);
        e = res.first; ree;
        combinedResults.append(res.second);
    }
#endif

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildUniqueMsInfoScanKeyVsScanPoints(
        const QString &msDataFilePath,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
        QMap<ScanNumber, ScanPoints > *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msDataFilePath); ree;

    e = msReaderParquet.collateTandemPrecursorTargetsDIA(diaTargetFrames); ree

    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderParquet.getMsScanInfos();
    for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {
        scanNumberVsScanTime->insert(it.key(), it.value().scanTime);
    }

    const int msLevel = 1;
    e = msReaderParquet.getScanPoints(msLevel, scanNumberVsScanTimeMS1); ree;

    ERR_RETURN
}


