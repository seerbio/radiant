//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretronProcessormatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "PSMsReader.h"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/serialization/vector.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using rTreeCoor = bg::model::point<int, 1, bg::cs::cartesian>;
using rTreePeakBox = bg::model::box<rTreeCoor>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreePeakBox, QString> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


Err MsFrameScoretronProcessormatic::init(
        const QString &frameScoreVecFilePath,
        const QString &frameExtractedScansFilePath,
        const PythiaParameters &pythiaParameters,
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        QPair<double, double> mzTargetStartStop
) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(frameScoreVecFilePath); ree;
    e = ErrorUtils::fileExists(frameExtractedScansFilePath); ree;

    m_params = pythiaParameters;
    m_scoreVectorsFilePath = frameScoreVecFilePath;
    m_frameExtractedScansFilePath = frameExtractedScansFilePath;

    e = MsFrame::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop,
            &m_msFrame
    ); ree;

    e = buildPepStrWModsVsExtractedScanRow(); ree;
    e = buildFrameIndexVsApexScorePeptideStringWithMods(); ree;
    e = buildPepStringWithModsVsMS2Ions(); ree;
    e = rescoreCandidatesWithFullPrediction(); ree;

    const int precision = 2;
    const int maxIters = 50;
    const double stopTolerance = 1e-8;
    e = m_deconvolvotron.init(
            precision,
            m_params.mzMaxDataStructure,
            maxIters,
            stopTolerance,
            m_params.pValThreshold
    ); ree;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::buildPepStrWModsVsExtractedScanRow() {

    ERR_INIT

    QVector<ExtractedScansReaderRow> extractedScansReaderRows;

    e = ParquetReader::read(
            m_frameExtractedScansFilePath,
            &extractedScansReaderRows
            ); ree;

    for (const ExtractedScansReaderRow &er : extractedScansReaderRows) {
        m_pepStrWModsVsExtractedScanRow.insert(er.peptideStringWithMods, er);
    }

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::buildFrameIndexVsApexScorePeptideStringWithMods() {

    ERR_INIT

    QVector<MsFrameScoreVectorReaderRow> msFrameScoreVectorReaderRows;

    e = ParquetReader::read(
            m_scoreVectorsFilePath,
            &msFrameScoreVectorReaderRows
            ); ree;

    for (const MsFrameScoreVectorReaderRow &fsrr : msFrameScoreVectorReaderRows) {
        m_frameIndexVsApexScorePeptideStringWithMods[fsrr.frameIndexMaxScore].push_back(fsrr.peptideStringWithMods);

        e = ErrorUtils::isNotEmpty(fsrr.scorePerFrameIndexOfTargetVec); ree;

        const double ogScore = fsrr.scorePerFrameIndexOfTargetVec.at(fsrr.frameIndexMaxScore);
        m_pepStrWModsVsOgScore.insert(fsrr.peptideStringWithMods, ogScore);
        m_pepStrWModsVsIsDecoy.insert(fsrr.peptideStringWithMods, fsrr.isDecoy);
        m_pepStrWModsVsCharge.insert(fsrr.peptideStringWithMods, fsrr.charge);
    }

    ERR_RETURN
}


namespace {

    Err getFrameIndexMs2Ions(
            double mzMin,
            double mzMax,
            const QList<ExtractedScansReaderRow> &extractedScansReaderRows,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanPreds
    ) {

        ERR_INIT

        for (const ExtractedScansReaderRow &row : extractedScansReaderRows) {

            QVector<MS2Ion> pepPred;
            e = ParallelUtils::zip(
                    row.mzSearched,
                    row.intensitySearched,
                    &pepPred
            ); ree;

            FragLibReader::filterMs2IonsByMz(
                    mzMin,
                    mzMax,
                    &pepPred
            );

            scanPreds->insert(row.peptideStringWithMods, pepPred);
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretronProcessormatic::buildPepStringWithModsVsMS2Ions() {

    ERR_INIT

    e = getFrameIndexMs2Ions(
            m_params.mzMinDataStructure,
            m_params.mzMaxDataStructure,
            m_pepStrWModsVsExtractedScanRow.values(),
            &m_pepStrWModsVsMS2Ions
    ); ree;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::rescoreCandidatesWithFullPrediction() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frameIndexVsApexScorePeptideStringWithMods); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsOgScore); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsMS2Ions); ree;

    for (auto it = m_frameIndexVsApexScorePeptideStringWithMods.begin(); it != m_frameIndexVsApexScorePeptideStringWithMods.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const ScanNumber scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndex);
        const ScanPoints &scanPoints = m_msFrame.getScanPointsByScanNumber(scanNumber);
        const QVector<PeptideStringWithMods> &PeptideStringWithModsVec = it.value();

        for (const PeptideStringWithMods &peptideStringWithMods : PeptideStringWithModsVec) {

            e = ErrorUtils::isTrue(m_pepStrWModsVsMS2Ions.contains(peptideStringWithMods)); ree;
            const QVector<MS2Ion> &ms2Ions = m_pepStrWModsVsMS2Ions.value(peptideStringWithMods);

            const ExtractPoints extractPoints = MsUtils::extractPointsFromPoints(
                    scanPoints,
                    ms2Ions,
                    m_params.ms2ExtractionWidthPPM
            );

            const QPair<QVector<double>, QVector<double>> intensitiesUnzipped
                    = ParallelUtils::unZip(extractPoints.intensityFoundVsSearched);

            const Eigen::VectorX<double> intensityFound = EigenUtils::convertQVectorToEigenVector(intensitiesUnzipped.first);
            const Eigen::VectorX<double> intensitySearched = EigenUtils::convertQVectorToEigenVector(intensitiesUnzipped.second);

            const QVector<QPointF> &mzFoundVsSearched = extractPoints.mzFoundVsSearched;

            const double cosineSim = EigenUtils::cosineSimilarity(intensityFound, intensitySearched);
            const double klDiv = EigenUtils::klDivergence(intensityFound, intensitySearched);

            const long ionsFound
                    = std::count_if(mzFoundVsSearched.begin(), mzFoundVsSearched.end(), [](const QPointF &p){return p.x() > 0;});

            const double fractionFound = static_cast<double>(ionsFound) / mzFoundVsSearched.size();

            const ReScore reScore = std::sqrt(MathUtils::pRound(
                    (cosineSim / klDiv)
                    * std::pow(ionsFound, 2)
                    * fractionFound
                    * std::pow(ionsFound, 3),
                    S_GLOBAL_SETTINGS.ROUNDING_PRECISION)
                            );

            ReScoreVals reScoreVals;
            reScoreVals.reScore = reScore;
            reScoreVals.klDiv = klDiv < 0 ? 1e3 : klDiv;
            reScoreVals.cosineSim = cosineSim;
            reScoreVals.fractionFound = fractionFound;
            reScoreVals.ionsFound = static_cast<int>(ionsFound);
            reScoreVals.score = m_pepStrWModsVsOgScore.value(peptideStringWithMods);

            m_pepStringModsVsRescore.insert(peptideStringWithMods, reScoreVals);
        }

    }

    ERR_RETURN
}

namespace {

    void sortDiscScoresDesc(
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>>> *frameIndexVsPeptideStringWithModsDiscScore
    ) {

        for (auto it = frameIndexVsPeptideStringWithModsDiscScore->begin();
             it != frameIndexVsPeptideStringWithModsDiscScore->end();
             it++
                ) {

            QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>> &scores = it.value();
            using Sort = QPair<PeptideStringWithMods, TandemDeconvolverResult>;
            std::sort(scores.rbegin(), scores.rend(), [](const Sort &l, const Sort &r){
                return l.second.discScore < r.second.discScore;
            });

        }

    }

}//namespace
Err MsFrameScoretronProcessormatic::processFrameScoreVectors(QVector<PSMsReaderRow> *psmReaderRows) {

    ERR_INIT

    psmReaderRows->clear();

    e = ErrorUtils::isTrue(m_msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsExtractedScanRow); ree;
    e = ErrorUtils::isNotEmpty(m_frameIndexVsApexScorePeptideStringWithMods); ree;

    for (auto it = m_frameIndexVsApexScorePeptideStringWithMods.begin();
         it != m_frameIndexVsApexScorePeptideStringWithMods.end();
         it++
        ) {

        const FrameIndex frameIndex = it.key();
        const QVector<PeptideStringWithMods> &peptideStringWithMods = it.value();

        if (peptideStringWithMods.isEmpty()) {
            continue;
        }

        e = processorLogic(peptideStringWithMods, frameIndex); ree;
    }

    sortDiscScoresDesc(&m_frameIndexVsPeptideDeconResult);

    e = compileScores(psmReaderRows); ree;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::processorLogic(
        const QVector<PeptideStringWithMods> &peptideStringWithMods,
        FrameIndex frameIndex
        ) {

    ERR_INIT

    const ScanNumber &scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndex);
    const ScanPoints &scanPoints = m_msFrame.getScanPointsByScanNumber(scanNumber);

    QMap<PeptideStringWithMods, QVector<MS2Ion>> peptideByExtractedPoints;
    e = collateMS2IonsByPepStrWithModsForFrameIndex(
            peptideStringWithMods,
            &peptideByExtractedPoints
            ); ree;

    e = calculateDiscriminateScoreForFrameIndex(
            peptideByExtractedPoints,
            scanPoints,
            frameIndex
    ); ree;

    e = buildFrameIndexVsFrameStats(); ree;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::collateMS2IonsByPepStrWithModsForFrameIndex(
        const QVector<PeptideStringWithMods> &peptideStringWithMods,
        QMap<PeptideStringWithMods, QVector<MS2Ion>> *peptideByExtractedPoints
                ) {

    ERR_INIT

    using SRT = QPair<PeptideStringWithMods, ReScore>;

    QVector<SRT> sortingPairVec;
    for (const PeptideStringWithMods &pswm : peptideStringWithMods) {
        e = ErrorUtils::isTrue(m_pepStringModsVsRescore.contains(pswm)); ree;
        const ReScore reScore = m_pepStringModsVsRescore.value(pswm).reScore;
        sortingPairVec.push_back({pswm, reScore});
    }

    const auto sortDescScore = [](const SRT &l, const SRT &r){
        return l.second < r.second;
    };

    std::sort(sortingPairVec.rbegin(), sortingPairVec.rend(), sortDescScore);
    sortingPairVec.resize(std::min(m_params.returnPSMTopN, sortingPairVec.size()));

    for (const SRT &pr : sortingPairVec) {
        const PeptideStringWithMods &pswm = pr.first;
        peptideByExtractedPoints->insert(pswm, m_pepStrWModsVsMS2Ions.value(pswm));
    }

    ERR_RETURN
}

namespace {

    Err buildUniqueMzVals(
            const QList<QVector<MS2Ion>> &ms2Ions,
            QVector<double> *ms2IonsUnique
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        const int precision = 4;

        QSet<int> mzHashedVals;

        for (const QVector<MS2Ion> &ions : ms2Ions) {

            for (const MS2Ion &ion : ions) {
                const int mzHashed = MathUtils::hashDecimal(ion.x(), precision);
                mzHashedVals.insert(mzHashed);
            }

            for (int mzHashed : mzHashedVals) {
                const auto mzUnhashed = MathUtils::unHashDecimal<double>(mzHashed, precision);
                ms2IonsUnique->push_back(mzUnhashed);
            }
        }

        ERR_RETURN
    }

    Err removeNegativeScanPredsFromFirstPassTandemDeconResults(
            const QMap<PeptideStringWithMods, TandemDeconvolverResult> &pepSeqVsWeight,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &scanPreds,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanPredsFiltered
            ) {

        ERR_INIT

        for (auto it = pepSeqVsWeight.begin(); it != pepSeqVsWeight.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const TandemDeconvolverResult &tandemDeconvolverResult = it.value();

            if (tandemDeconvolverResult.discScore < 0) {
                continue;
            }

            e = ErrorUtils::isTrue(scanPreds.contains(peptideStringWithMods)); ree;

            scanPredsFiltered->insert(peptideStringWithMods, scanPreds.value(peptideStringWithMods));
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretronProcessormatic::calculateDiscriminateScoreForFrameIndex(
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &peptideByExtractedPoints,
        const ScanPoints &scanPoints,
        FrameIndex frameIndex
) {

    ERR_INIT

    QVector<double> mzValsUnique;
    e = buildUniqueMzVals(
            peptideByExtractedPoints.values(),
            &mzValsUnique
            ); ree;

    ScanPoints extractedScanPoints = MsUtils::extractPointsFromPoints(
            scanPoints,
            mzValsUnique,
            m_params.ms2ExtractionWidthPPM
            );

    if (extractedScanPoints.isEmpty()) {
        extractedScanPoints.push_back({0, 1});
    }

    QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsTandemDeconvolverResult;
    e = m_deconvolvotron.deconvolveTandemSpectra(
            extractedScanPoints,
            peptideByExtractedPoints,
            &pepSeqVsTandemDeconvolverResult
    ); ree;


//#define REMOVE_NEG_DISC_SCORES
#ifdef REMOVE_NEG_DISC_SCORES
    QMap<PeptideStringWithMods, QVector<MS2Ion>> scanPredsFiltered;
    e = removeNegativeScanPredsFromFirstPassTandemDeconResults(
            pepSeqVsTandemDeconvolverResult,
            peptideByExtractedPoints,
            &scanPredsFiltered
            ); ree;

    e = m_deconvolvotron.deconvolveTandemSpectra(
            extractedScanPoints,
            scanPredsFiltered,
            &pepSeqVsTandemDeconvolverResult
    ); ree;
#endif

    for (auto itt = pepSeqVsTandemDeconvolverResult.begin(); itt != pepSeqVsTandemDeconvolverResult.end(); itt++) {

        const PeptideStringWithMods &peptideStringWithMods = itt.key();
        const TandemDeconvolverResult discriminateScore = itt.value();

        m_frameIndexVsPeptideDeconResult[frameIndex].push_back({peptideStringWithMods, discriminateScore});
    }

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::buildFrameIndexVsFrameStats() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frameIndexVsPeptideDeconResult); ree;
    e = ErrorUtils::isNotEmpty(m_frameIndexVsApexScorePeptideStringWithMods); ree;
    e = ErrorUtils::isNotEmpty(m_pepStringModsVsRescore); ree;

    for (auto it = m_frameIndexVsApexScorePeptideStringWithMods.begin();
        it != m_frameIndexVsApexScorePeptideStringWithMods.end();
        it++
        ) {

        const FrameIndex frameIndex = it.key();
        const QVector<PeptideStringWithMods> &peptideStringWithMods = it.value();
        
        QVector<double> reScores;
        reScores.reserve(peptideStringWithMods.size());
        
        for (const PeptideStringWithMods &pep : peptideStringWithMods) {

            e = ErrorUtils::isTrue(m_pepStringModsVsRescore.contains(pep)); ree;

            const ReScoreVals &reScoreVals = m_pepStringModsVsRescore.value(pep);
            reScores.push_back(reScoreVals.reScore);
        }

        QVector<double> discScores;
        discScores.reserve(peptideStringWithMods.size());
        
        const QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>> &pr 
                = m_frameIndexVsPeptideDeconResult.value(frameIndex);
        
        std::transform(
                pr.begin(), 
                pr.end(), 
                std::back_inserter(discScores),
                [](const QPair<PeptideStringWithMods, TandemDeconvolverResult> &r){return r.second.discScore;}
                );
        
        const auto discScoreMinMax = std::minmax_element(discScores.begin(), discScores.end());
        const auto reScoreMinMax = std::minmax_element(reScores.begin(), reScores.end());
        
        FrameStats fs;
        fs.discScoreMin = *discScoreMinMax.first;
        fs.discScoreMax = *discScoreMinMax.second;
        fs.discScoreMean = MathUtils::mean(discScores);
        fs.discScoreMedian = MathUtils::median(discScores);
        fs.discScoreStDev = MathUtils::stDev(discScores);

        fs.scoreMin = *reScoreMinMax.first;
        fs.scoreMax = *reScoreMinMax.second;
        fs.scoreMean = MathUtils::mean(reScores);
        fs.scoreMedian = MathUtils::median(reScores);
        fs.scoreStDev = MathUtils::stDev(reScores);

        fs.frameCandidateCount = peptideStringWithMods.size();

        m_frameIndexVsPeptideFrameStats.insert(frameIndex, fs);
    }
    
    ERR_RETURN
}

namespace {

    //TODO make this also viable for n-terminal cleave points
    int calculateMissedCleavages(
            const QString &peptideStringWithMods,
            const QStringList &cTermCleavePoints
            ) {

        int cleavePointCount = 0;

        for (const QString &cleavePoint : cTermCleavePoints) {
            cleavePointCount += peptideStringWithMods.count(cleavePoint);
        }

        cleavePointCount = cTermCleavePoints.contains(peptideStringWithMods.back())
                ? cleavePointCount - 1
                : cleavePointCount;

        return std::max(0, cleavePointCount);
    }

}
Err MsFrameScoretronProcessormatic::compileScores(QVector<PSMsReaderRow> *psmReaderRows) {
    
    ERR_INIT
    
    e = ErrorUtils::isNotEmpty(m_frameIndexVsPeptideDeconResult); ree;
    e = ErrorUtils::isNotEmpty(m_pepStringModsVsRescore); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsOgScore); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsIsDecoy); ree;
    e = ErrorUtils::isNotEmpty(m_pepStrWModsVsCharge); ree;
    e = ErrorUtils::isAboveThreshold(m_msFrame.scanCount(), 0, ErrorUtilsParam::ExcludeThreshold); ree;

    using SRT = QPair<PeptideStringWithMods, TandemDeconvolverResult>;
    const auto sortDiscScoresDescLogic
            = [](const SRT &l, const SRT &r){return l.second.discScore < r.second.discScore;};

    const auto sortReScoresDescLogic
            = [](const PSMsReaderRow &l, const PSMsReaderRow &r){return l.rescore < r.rescore;};
    
    for (auto it = m_frameIndexVsPeptideDeconResult.begin(); 
            it != m_frameIndexVsPeptideDeconResult.end(); 
            it++) {
        
        const FrameIndex frameIndex = it.key();
        QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>> prs = it.value();

        std::sort(prs.rbegin(), prs.rend(), sortDiscScoresDescLogic);

        const FrameStats &frameStats = m_frameIndexVsPeptideFrameStats.value(frameIndex);
        QVector<PSMsReaderRow> psmReaderRowsFrameIndex;

        int discScoreCounter = 0;
        for (const QPair<PeptideStringWithMods, TandemDeconvolverResult> &pr : prs) {
            
            const PeptideStringWithMods &peptideStringWithMods = pr.first;
            const TandemDeconvolverResult &deconRes = pr.second;
            const ReScoreVals &reScoreVals = m_pepStringModsVsRescore.value(peptideStringWithMods);
            
            PSMsReaderRow row;
            row.frameIndex = frameIndex;
            row.scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndex);
            row.charge = m_pepStrWModsVsCharge.value(peptideStringWithMods);
            row.uniqueMsInfoScanKey = m_msFrame.uniqueMsInfoScanKey();

            row.peptideStringWithMods = peptideStringWithMods;

            row.score = reScoreVals.score;

            row.discScore = deconRes.discScore;
            row.frameRankDiscScore = discScoreCounter++;

            row.isDecoy = m_pepStrWModsVsIsDecoy.value(peptideStringWithMods);
            row.rescore = reScoreVals.reScore;

            row.pVal = deconRes.pVal;
            row.tTest = deconRes.tTestVal;
            row.frameFStat = deconRes.frameFStat;
            row.pValFrameFtest = deconRes.pValFrameFtest;
            row.frameError = deconRes.frameError;

            row.scoreMean = frameStats.scoreMean;
            row.scoreMin = frameStats.scoreMin;
            row.scoreMax = frameStats.scoreMax;
            row.scoreMedian = frameStats.scoreMedian;
            row.scoreStDev = frameStats.scoreStDev;

            row.discScoreMean = frameStats.discScoreMean;
            row.discScoreMin = frameStats.discScoreMin;
            row.discScoreMax = frameStats.discScoreMax;
            row.discScoreMedian = frameStats.discScoreMedian;
            row.discScoreStDev = frameStats.discScoreStDev;

            row.frameCandidateCount = frameStats.frameCandidateCount;
            
            row.cosineSim = reScoreVals.cosineSim;
            row.fractionFound = reScoreVals.fractionFound;
            row.klDiv = reScoreVals.klDiv;
            row.ionsFound = reScoreVals.ionsFound;

            row.peptideSize = row.peptideStringWithMods.size();

            row.missedCleavages = calculateMissedCleavages(
                    row.peptideStringWithMods,
                    m_params.cTermCleavePoints
                    );

            psmReaderRowsFrameIndex.push_back(row);
        }

        std::sort(psmReaderRowsFrameIndex.rbegin(), psmReaderRowsFrameIndex.rend(), sortReScoresDescLogic);

        int reScoreCounter = 0;
        for(PSMsReaderRow &psmRow : psmReaderRowsFrameIndex) {
            psmRow.frameRankScore = reScoreCounter++;
            psmReaderRows->push_back(psmRow);
        }

    }
    
    ERR_RETURN
}

