//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FrameExtractsReader.h"
#include "MsFrameScoretronProcessormatic.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>


namespace {

    Err setFeatureFinderParams(
            const PythiaParameters &pythiaParameters,
            FeatureFinderParameters *featureFinderParameters
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        featureFinderParameters->sigma = pythiaParameters.sigma;
        featureFinderParameters->filterLength = pythiaParameters.filterLength;
        featureFinderParameters->smoothCount = pythiaParameters.smoothCount;
        featureFinderParameters->signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        featureFinderParameters->tolerancePPM = pythiaParameters.ms2ExtractionWidthPPM;
        featureFinderParameters->skipScanCount = pythiaParameters.skipScanCount;
        featureFinderParameters->minScanCount = pythiaParameters.minScanCount;

        featureFinderParameters->printParams();

        e = ErrorUtils::isTrue(featureFinderParameters->isValid()); ree

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const QString &fragLibBackgroundFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    params.print();

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isAboveThreshold(
            mzTargetStartStop.second,
            mzTargetStartStop.first,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;

    m_params = params;
    m_msDataFilePath = msDataFilePath;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_mzTargetStartStop = mzTargetStartStop;

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibFilePath,
            m_params.chargeStateMin,
            m_params.chargeStateMax,
            m_mzTargetStartStop.first,
            m_mzTargetStartStop.second,
            &m_fragPreds,
            &m_fragPredsIsDecoy
    ); ree

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibBackgroundFilePath,
            m_params.chargeStateMin,
            m_params.chargeStateMax,
            m_mzTargetStartStop.first,
            m_mzTargetStartStop.second,
            &m_fragPredsBackground,
            &m_fragPredsBackgroundIsDecoy
    ); ree

    e = ErrorUtils::isNotEmpty(m_fragPreds); ree
    e = m_fragLibIonRTree.init(m_fragPreds); ree

    e = ErrorUtils::isNotEmpty(m_fragPredsBackground); ree
    e = m_fragLibIonRTreeBackground.init(m_fragPredsBackground); ree

    QMap<MzHashed, FrequencyPercent> mzHashVsFreqPctBackground;
    e = m_fragLibIonRTreeBackground.buildMzHashedVsFragLibIonFrequencePercentages(
            m_params.ms2ExtractionWidthPPM,
            &mzHashVsFreqPctBackground
            ); ree

    e = m_fragLibIonRTree.addFrequencyPercentagesToFragLibIons(mzHashVsFreqPctBackground); ree
    e = m_fragLibIonRTreeBackground.addFrequencyPercentagesToFragLibIons(mzHashVsFreqPctBackground); ree

    e = initMsFrame(
        msDataFilePath,
        uniqueMsInfoScanKey,
        mzTargetStartStop
        ); ree

    qDebug() << "TargetKey" << uniqueMsInfoScanKey;
    qDebug() << "Scan Count" << m_msFrame.scanCount();
    qDebug() << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}

Err MsFrameScoretron::initMsFrame(
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    e = MsFrame::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop,
            &m_msFrame
    ); ree;

    e = m_msFrame.filterFrameByMz(
            m_params.mzMinDataStructure,
            m_params.mzMaxDataStructure
    ); ree;

//    e = m_msFrame.smoothFrame(
//            m_params.filterLength,
//            m_params.sigma,
//            m_params.smoothCount,
//            m_params.mzMaxDataStructure
//    ); ree

    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
    );ree;




    ERR_RETURN
}

namespace {

    Err buildMsFrameApexTurboXIC(
            const PythiaParameters &pythiaParameters,
            const QMap<FrameIndex, ScanPoints> &frameIndexVsScanPoints,
            TurboXIC *turboXic
            ) {

        ERR_INIT

        FeatureFinderParameters featureFinderParameters;
        e = setFeatureFinderParams(pythiaParameters, &featureFinderParameters); ree

        FeatureFinderHillBuilder featureFinderHillBuilder;
        e = featureFinderHillBuilder.init(featureFinderParameters); ree

        e = featureFinderHillBuilder.buildHills(frameIndexVsScanPoints); ree

        QVector<FeatureFinderHill> featureFinderHills;
        e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills); ree

        //TODO consider making these settable in parameters
        const int smoothCount = 2;
        const int filterLength = 5;
        const int order = 1;
        const int derivative = 0;
        const int rate = 1;

        QMap<FrameIndex, ScanPoints> frameApexes;

        for (const FeatureFinderHill &ffh : featureFinderHills) {

            Eigen::VectorX<double> intensitiesVec = EigenUtils::convertQVectorToEigenVector(ffh.intensities());

            for (int i = 0; i < smoothCount; i++) {
                e = EigenKernelUtils::savitskyGolaySmooth(
                        filterLength,
                        order,
                        derivative,
                        rate,
                        &intensitiesVec
                        ); ree
            }

            const QVector<FrameIndex> &frameIndexes = ffh.scanNumbers();
            const double mzMean = ffh.mzMean();
            const QVector<double> &originalIntensities = ffh.intensities();

            const QMap<int, double> apexes = EigenUtils::apexes(intensitiesVec, 1);

            for (auto it = apexes.begin(); it != apexes.end(); it++) {

                const int apexIndex = it.key();
                const double _smoothedIntensityApexValue = it.value();
                const int frameIndex = frameIndexes.at(apexIndex);
                const double originalIntensity = originalIntensities.at(apexIndex);

                frameApexes[frameIndex].push_back({mzMean, originalIntensity});
            }

        }

        turboXic->init(frameApexes); ree

        ERR_RETURN

    }

    Err buildApexSpectra(
            const TurboXIC &frameApexesTurboXIC,
            QMap<FrameIndex, ScanPoints> *apexScanPoints
            ) {

        ERR_INIT

        double frameIndexMin;
        double frameIndexMax;
        double mzMin;
        double mzMax;

        e = frameApexesTurboXIC.getRTreeLimits(
                &frameIndexMin,
                &frameIndexMax,
                &mzMin,
                &mzMax
                ); ree

        for (
            auto frameIndex = static_cast<int>(frameIndexMin);
            frameIndex <= static_cast<int>(frameIndexMax);
            ++frameIndex
            ) {

            //TODO figure out if this is needed
            const int frameIndexLo = std::max(0, frameIndex - 1);
            const int frameIndexHi = std::min(static_cast<int>(frameIndexMax), frameIndex + 1);

            const ScanPoints frameIndexScanPoints = frameApexesTurboXIC.extractSpectrum(
                    mzMin,
                    mzMax,
                    frameIndex,
                    frameIndex
                    );

            apexScanPoints->insert(frameIndex, frameIndexScanPoints);
        }

        ERR_RETURN
    }


}//namespace
Err MsFrameScoretron::scoreFrameCandidates(QMap<FrameIndex , QVector<ScoredCandidate>> *frameIndexVsScoredCandidates) {

    ERR_INIT

    TurboXIC frameApexesTurboXIC;
    e = buildMsFrameApexTurboXIC(
            m_params,
            m_msFrame.frameIndexVsScanPoints(),
            &frameApexesTurboXIC
            ); ree

    QMap<FrameIndex, ScanPoints> apexScanPoints;
    e = buildApexSpectra(
            frameApexesTurboXIC,
            &apexScanPoints
            ); ree

//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    for(auto it = apexScanPoints.begin(); it != apexScanPoints.end(); it++) {
        const ScanNumber sn = m_msFrame.scanNumberFromFrameIndex(it.key());
        scanNumberVsScanPoints.insert(sn, it.value());
    }

    e = MsFrame::writeFrameScans(scanNumberVsScanPoints, "frame.prq"); ree
#endif

    e = iterateApexScanPoints(
            apexScanPoints,
            frameIndexVsScoredCandidates
            ); ree


    ERR_RETURN
}

namespace {

    QMap<PeptideId, QVector<FragLibIon>> filterPeptideIdVsFragLibIonsForFrameIndex(
            const QMap<PeptideId, QVector<FragLibIon>> &peptideIdVsFragLibIonsForFrameIndex,
            int fragLibIonCountForCandidateMin
            ) {

        QMap<PeptideId, QVector<FragLibIon>> peptideIdVsFragLibIonsForFrameIndexFiltered;

        for (
            auto it = peptideIdVsFragLibIonsForFrameIndex.begin();
            it != peptideIdVsFragLibIonsForFrameIndex.end();
            it++
            ) {

            const PeptideId peptideId = it.key();
            const QVector<FragLibIon> &flis = it.value();

            if (flis.size() < fragLibIonCountForCandidateMin) {
                continue;
            }

            peptideIdVsFragLibIonsForFrameIndexFiltered.insert(peptideId, flis);
        }

        return peptideIdVsFragLibIonsForFrameIndexFiltered;
    }

    QVector<FragLibIon> removeIsotopeFragLibsWithoutMatchingMonoIsotope(const QVector<FragLibIon> &foundFragLibIons) {

        QVector<FragLibIon> foundFragLibIonsSortedByIsIsotope = foundFragLibIons;

        const auto sortLogicIsIso = [](const FragLibIon &l, const FragLibIon &r){return l.isIsotope < r.isIsotope;};
        std::sort(foundFragLibIonsSortedByIsIsotope.begin(), foundFragLibIonsSortedByIsIsotope.end(), sortLogicIsIso);

        QHash<QString, bool> monoIsotopeIndexes;
        QVector<FragLibIon> filteredIsotopes;
        for (const FragLibIon &fli : foundFragLibIonsSortedByIsIsotope) {

            const QString indexHash = fli.ionType + QString::number(fli.ionIndex);

            if (!fli.isIsotope) {
                filteredIsotopes.push_back(fli);
                monoIsotopeIndexes.insert(indexHash, true);
                continue;
            }

            if (monoIsotopeIndexes.value(indexHash)) {
                filteredIsotopes.push_back(fli);
            }
        }

        const auto sortLogicMzAsc = [](const FragLibIon &l, const FragLibIon &r){return l.mz < r.mz;};
        std::sort(filteredIsotopes.begin(), filteredIsotopes.end(), sortLogicMzAsc);

        return filteredIsotopes;
    }

    double sumFragLibIonFrequencyPercent(const QVector<FragLibIon> &fragLibIons) {
        const auto sumLogic = [](double sum, const FragLibIon &fli){
            return sum + fli.frequencyPercent;
        };
        return std::accumulate(fragLibIons.begin(), fragLibIons.end(), 0.0, sumLogic);
    }

    QPair<double, double> calcKLDivergenceCosineSim(
            const QVector<FragLibIon> &tandemPredictionFragLibIons,
            const QVector<FragLibIon> &foundFragLibIonsFiltered
            ) {

        QHash<QString, QPair<double, double>> fragLibIonAlignment;

        for (const FragLibIon &fli : tandemPredictionFragLibIons) {

            const QString key = fli.ionType + QString::number(fli.ionIndex) + QString::number(std::round(fli.mz * 10));
            fragLibIonAlignment[key] = {fli.intensity, 0.0};
        }

        for (const FragLibIon &fli : foundFragLibIonsFiltered) {

            const QString key = fli.ionType + QString::number(fli.ionIndex) + QString::number(std::round(fli.mz * 10));
            fragLibIonAlignment[key].second = fli.intensitySearched;
        }

        Eigen::VectorX<double> predictionVec(fragLibIonAlignment.size());
        Eigen::VectorX<double> foundVec(fragLibIonAlignment.size());

        int index = 0;
        for (const QPair<double, double> &pr : fragLibIonAlignment) {
            predictionVec.coeffRef(index) = pr.first;
            foundVec.coeffRef(index) = pr.second;
            index++;
        }

        const double klDiv = EigenUtils::klDivergence(predictionVec, foundVec);
        const double cosSim = EigenUtils::cosineSimilarity(predictionVec, foundVec);

        return {klDiv, cosSim};
    }

}//namespace
Err MsFrameScoretron::iterateApexScanPoints(
        const QMap<FrameIndex, ScanPoints> &apexScanPoints,
        QMap<FrameIndex , QVector<ScoredCandidate>> *frameIndexVsScoredCandidates
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(apexScanPoints); ree

    for (auto it = apexScanPoints.begin(); it != apexScanPoints.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const ScanPoints &scanPoints = it.value();

        QMap<PeptideId, QVector<FragLibIon>> peptideIdVsFragLibIonsCandidatesForFrameIndex;
        e = extractFragLibIonsForScanPoints(
                scanPoints,
                &m_fragLibIonRTree,
                &peptideIdVsFragLibIonsCandidatesForFrameIndex
                ); ree

        QVector<ScoredCandidate> scoredCandidatesTarget;
        e = extractScores(
                peptideIdVsFragLibIonsCandidatesForFrameIndex,
                &scoredCandidatesTarget
        ); ree

        frameIndexVsScoredCandidates->insert(frameIndex, scoredCandidatesTarget);

        qDebug() << "processed frame index:" << frameIndex;
    }

    ERR_RETURN
}

Err MsFrameScoretron::extractFragLibIonsForScanPoints(
        const ScanPoints &scanPoints,
        FragLibIonRTree *fragLibIonRTree,
        QMap<PeptideId, QVector<FragLibIon>> *peptideIdVsFragLibIonsForFrameIndexOutput
        ) const {

    ERR_INIT

    QMap<PeptideId, QVector<FragLibIon>> peptideIdVsFragLibIonsForFrameIndex;

    for (const ScanPoint &sp : scanPoints) {

        if (sp.x() < m_params.mzMinDataStructure) {
            continue;
        }

        const double mzTol = MathUtils::calculatePPM(sp.x(), m_params.ms2ExtractionWidthPPM);
        const double mzMin = sp.x() - mzTol;
        const double mzMax = sp.x() + mzTol;

        QVector<FragLibIon> foundFragLibIonsForMz;

        // TODO when iRT is incorporated, use other overloaded function that takes iRT min/max as args
        e = fragLibIonRTree->getFragLibIons(
                mzMin,
                mzMax,
                &foundFragLibIonsForMz
        ); ree

        for (FragLibIon &fli : foundFragLibIonsForMz) {
            fli.mzSearched = sp.x();
            fli.intensitySearched = sp.y();
            peptideIdVsFragLibIonsForFrameIndex[fli.peptideId].push_back(fli);
        }
    }

    *peptideIdVsFragLibIonsForFrameIndexOutput = filterPeptideIdVsFragLibIonsForFrameIndex(
            peptideIdVsFragLibIonsForFrameIndex,
            m_params.minFoundMzPeaks
    );

    ERR_RETURN
}


Err MsFrameScoretron::extractScores(
        const QMap<PeptideId, QVector<FragLibIon>> &peptideIdVsFragLibIonsForFrameIndex,
        QVector<ScoredCandidate> *scoredCandidates
) {

    ERR_INIT

    for (auto it = peptideIdVsFragLibIonsForFrameIndex.begin(); it != peptideIdVsFragLibIonsForFrameIndex.end(); it++) {

        const PeptideId peptideId = it.key();
        const QVector<FragLibIon> &foundFragLibIons = it.value();

        QVector<FragLibIon> tandemPredictionFragLibIons;
        e = m_fragLibIonRTree.getFragLibIonsByPeptideId(
                peptideId,
                &tandemPredictionFragLibIons
        ); ree

        const double bestPossibleFreqPctScore = sumFragLibIonFrequencyPercent(tandemPredictionFragLibIons);

        const QVector<FragLibIon> foundFragLibIonsFiltered
                = removeIsotopeFragLibsWithoutMatchingMonoIsotope(foundFragLibIons);

        const double foundFreqPctScore = sumFragLibIonFrequencyPercent(foundFragLibIonsFiltered);

        const QPair<double, double> klDivergenceVsCosineSim = calcKLDivergenceCosineSim(
                tandemPredictionFragLibIons,
                foundFragLibIonsFiltered
        );

        ScoredCandidate scoredCandidate;
        scoredCandidate.frequencyPercentSum = foundFreqPctScore;
        scoredCandidate.frequencyPercentSumBestPossible = bestPossibleFreqPctScore;
        scoredCandidate.klDivergence = klDivergenceVsCosineSim.first;
        scoredCandidate.cosineSim = klDivergenceVsCosineSim.second;
        e = m_fragLibIonRTree.getPeptideSequenceWithMods(peptideId, &scoredCandidate.peptideStringWithMods); ree
        scoredCandidate.isDecoy = m_fragPredsIsDecoy.value(scoredCandidate.peptideStringWithMods);

        scoredCandidates->push_back(scoredCandidate);
    }

    ERR_RETURN
}