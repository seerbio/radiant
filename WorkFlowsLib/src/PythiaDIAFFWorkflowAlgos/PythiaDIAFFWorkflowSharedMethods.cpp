//
// Created by andrewnichols on 9/28/24.
//

#include "PythiaDIAFFWorkflowSharedMethods.h"

#include "DiscriminantScoretron.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FDRCLassifierNeuralNet.h"
#include "MS2Ion.h"
#include "MsFrame.h"
#include "QValueSettertron.h"
#include "TurboXIC.h"

#include <algorithm>
#include <cmath>
#include <QSet>


Err PythiaDIAFFWorkflowSharedMethods::buildUniqueMsScanInfosForProcessing(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        int numberOfUniqueScanInfosForPurpose,
        int offset,
        QVector<MsScanInfo> *uniqueMsScanInfosPurpose
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
    e = ErrorUtils::isAboveThreshold(numberOfUniqueScanInfosForPurpose, 8, ErrorUtilsParam::IncludeThreshold); ree;

    const double halfSizeScanInfos = uniqueMsScanInfos.size() / 2.0;
    const int incrementSize = std::max(1, static_cast<int>(halfSizeScanInfos / numberOfUniqueScanInfosForPurpose));

    for (int i = 0; i < uniqueMsScanInfos.size(); i += incrementSize) {

        const int distanceFromCenterRight
                = std::min(static_cast<int>(halfSizeScanInfos) + i + offset, uniqueMsScanInfos.size() - 1);

        const int distanceFromCenterLeft = std::max(static_cast<int>(halfSizeScanInfos) - i - offset, 0);

        if (i == 0) {
            uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
            continue;
        }

        uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
        uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterLeft));

        if (uniqueMsScanInfosPurpose->size() >= numberOfUniqueScanInfosForPurpose) {
            break;
        }

    }

    uniqueMsScanInfosPurpose->resize(std::min(
        numberOfUniqueScanInfosForPurpose,
        uniqueMsScanInfosPurpose->size()
        ));

    ERR_RETURN
}

namespace {

    constexpr float TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0 = 0.03f;
    constexpr float TIMS_CALIBRATED_IM_FILTER_PAD_MIN_ONE_OVER_K0 = 0.01f;
    constexpr float TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0 = 0.15f;
    constexpr int TIMS_MAIN_EVIDENCE_TOP_N_MS2_IONS = 4;

    struct TimsEvidenceCandidate {
        TargetDecoyCandidatePair *targetDecoyPair = nullptr;
        double evidenceScore = 0.0;
    };

    float ionMobilityCenterForFilter(
        const TargetDecoyCandidatePair *candidate,
        const MsCalibratomatic *msCalibratomatic
        ) {

        if (candidate == nullptr) {
            return -1.0f;
        }

        const float libraryIonMobility = candidate->iIM();
        if (libraryIonMobility <= 0.0f
            || msCalibratomatic == nullptr
            || !msCalibratomatic->isInitIM()) {
            return libraryIonMobility;
        }

        float predictedIonMobility = 0.0f;
        if (msCalibratomatic->predictIonMobility(
            libraryIonMobility,
            &predictedIonMobility
            ) != eNoError) {
            return libraryIonMobility;
        }

        return predictedIonMobility > 0.0f
             ? predictedIonMobility
             : libraryIonMobility;
    }

    float ionMobilityFilterPad(
        const PythiaParameters &pythiaParameters,
        const MsCalibratomatic *msCalibratomatic
        ) {

        if (msCalibratomatic == nullptr || !msCalibratomatic->isInitIM()) {
            return TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        }

        const float calibratedPad = msCalibratomatic->ionMobilityStDev(
            pythiaParameters.scanTimeWindowStDevs
            );
        if (calibratedPad <= 0.0f) {
            return TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        }

        return std::max(
            TIMS_CALIBRATED_IM_FILTER_PAD_MIN_ONE_OVER_K0,
            std::min(calibratedPad, TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0)
            );
    }

    bool isLibraryIonMobilityInAcquisitionWindow(
        const TargetDecoyCandidatePair *candidate,
        const PythiaParameters &pythiaParameters,
        const MsCalibratomatic *msCalibratomatic,
        const MsScanInfo &msScanInfo
        ) {

        if (candidate == nullptr) {
            return false;
        }

        const float ionMobilityCenter = ionMobilityCenterForFilter(
            candidate,
            msCalibratomatic
            );
        if (ionMobilityCenter <= 0.0f || msScanInfo.ionMobilityDriftTime <= 0.0f) {
            return true;
        }

        float windowLower = msScanInfo.ionMobilityWindowLower;
        float windowUpper = msScanInfo.ionMobilityWindowUpper;
        if (windowLower <= 0.0f || windowUpper <= 0.0f || windowUpper < windowLower) {
            windowLower = msScanInfo.ionMobilityDriftTime - TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
            windowUpper = msScanInfo.ionMobilityDriftTime + TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
        }

        const float windowPad = ionMobilityFilterPad(pythiaParameters, msCalibratomatic);
        windowLower -= windowPad;
        windowUpper += windowPad;
        return windowLower <= ionMobilityCenter && ionMobilityCenter <= windowUpper;
    }

    struct TurboXICLoadInput {
        QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        MzTargetKey mzTargetKey;
    };

    std::tuple<Err, MzTargetKey, TurboXIC*> buildTurboXICLogic(const TurboXICLoadInput &turboXicLoadInput) {

        ERR_INIT

        MsFrame msFrame;
        e = msFrame.init(
            turboXicLoadInput.scanNumberVsScanPoints,
            turboXicLoadInput.scanNumberVsScanTime
            ); rtee;

        auto turboXic = new TurboXIC();
        e = turboXic->init(msFrame.frameIndexVsScanPoints()); rtee;

        return {e, turboXicLoadInput.mzTargetKey, turboXic};
    }

    double maxIntensityInFrameWindow(
        const TurboXIC &turboXic,
        const MS2Ion &ms2Ion,
        float ppmTolerance,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax
        ) {

        if (ms2Ion.mz <= 0.0f || ppmTolerance <= 0.0f) {
            return 0.0;
        }

        if (frameIndexMax < frameIndexMin) {
            std::swap(frameIndexMin, frameIndexMax);
        }

        const float massTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTolerance);
        const XICPoints xicPoints = turboXic.extractPointsXIC(ms2Ion.mz - massTol, ms2Ion.mz + massTol);

        double maxIntensity = 0.0;
        for (const XICPoint &xicPoint : xicPoints) {
            if (xicPoint.scanNumber < frameIndexMin || xicPoint.scanNumber > frameIndexMax) {
                continue;
            }
            if (xicPoint.intensity <= maxIntensity) {
                continue;
            }
            maxIntensity = static_cast<double>(xicPoint.intensity);
        }

        return maxIntensity;
    }

    bool predictedFrameWindow(
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const MsFrame &msFrameMzTarget,
        const TargetDecoyCandidatePair *targetDecoyPair,
        FrameIndex *frameIndexMin,
        FrameIndex *frameIndexMax
        ) {

        if (targetDecoyPair == nullptr
            || !msCalibratomatic.isInitRT()
            || !msFrameMzTarget.isValid()) {
            return false;
        }

        float predictedScanTime = -1.0f;
        Err e = msCalibratomatic.predictScanTime(targetDecoyPair->iRt(false), &predictedScanTime);
        if (e != eNoError) {
            return false;
        }

        const float scanTimeWindow
            = msCalibratomatic.scanTimeStDev(static_cast<float>(pythiaParameters.scanTimeWindowStDevs));

        e = msFrameMzTarget.frameIndexFromScanTime(predictedScanTime - scanTimeWindow, frameIndexMin);
        if (e != eNoError) {
            return false;
        }

        e = msFrameMzTarget.frameIndexFromScanTime(predictedScanTime + scanTimeWindow, frameIndexMax);
        if (e != eNoError) {
            return false;
        }

        if (*frameIndexMax < *frameIndexMin) {
            std::swap(*frameIndexMin, *frameIndexMax);
        }

        return true;
    }

    double calculateCheapTimsEvidenceScore(
        const TurboXIC &turboXic,
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const MsFrame &msFrameMzTarget,
        const TargetDecoyCandidatePair *targetDecoyPair
        ) {

        if (targetDecoyPair == nullptr) {
            return 0.0;
        }

        FrameIndex frameIndexMin = 0;
        FrameIndex frameIndexMax = 0;
        if (!predictedFrameWindow(
                msCalibratomatic,
                pythiaParameters,
                msFrameMzTarget,
                targetDecoyPair,
                &frameIndexMin,
                &frameIndexMax
                )) {
            return 0.0;
        }

        QVector<MS2Ion> ms2Ions = targetDecoyPair->ms2IonsTarget();
        if (ms2Ions.isEmpty()) {
            return 0.0;
        }

        MS2Ion::sortMS2IonsIntensityDesc(&ms2Ions);

        const int ionCount = std::min(TIMS_MAIN_EVIDENCE_TOP_N_MS2_IONS, ms2Ions.size());
        double weightedSignal = 0.0;
        int matchedIonCount = 0;
        for (int ionIndex = 0; ionIndex < ionCount; ++ionIndex) {
            const MS2Ion &ms2Ion = ms2Ions.at(ionIndex);
            const double maxIntensity = maxIntensityInFrameWindow(
                turboXic,
                ms2Ion,
                static_cast<float>(pythiaParameters.ms2ExtractionWidthPPM),
                frameIndexMin,
                frameIndexMax
                );
            if (maxIntensity <= 0.0) {
                continue;
            }

            const double libraryWeight = std::sqrt(std::max(1.0, static_cast<double>(ms2Ion.intensity)));
            const double signalWeight = std::log1p(maxIntensity) * libraryWeight;
            weightedSignal += signalWeight;
            ++matchedIonCount;
        }

        if (matchedIonCount == 0) {
            return 0.0;
        }

        return weightedSignal * static_cast<double>(matchedIonCount);
    }
}//namespace
Err PythiaDIAFFWorkflowSharedMethods::buildMzTargetKeyVsTurboXicPntrs(
    const QVector<MsScanInfo>& uniqueMsScanInfos,
    const QMap<ScanNumber, ScanTime>& scanNumberVsScanTime,
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
    QMap<MzTargetKey, TurboXIC*>* mzTargetKeyVsTurboXicPntr
    ) {

    ERR_INIT

            e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

    QVector<TurboXICLoadInput> turboXicLoadInputs;
    for (const MsScanInfo &msScanInfo : uniqueMsScanInfos) {

        const MzTargetKey mzTargetKey = msScanInfo.targetKey();
        e = ErrorUtils::contains(mzTargetKey, *diaTargetFrames); ree;

        TurboXICLoadInput turboXicLoadInput;
        turboXicLoadInput.scanNumberVsScanPoints = diaTargetFrames->value(mzTargetKey);
        turboXicLoadInput.scanNumberVsScanTime = scanNumberVsScanTime;
        turboXicLoadInput.mzTargetKey = mzTargetKey;

        turboXicLoadInputs.push_back(turboXicLoadInput);
    }

    QFuture<std::tuple<Err, MzTargetKey, TurboXIC*>> futures = QtConcurrent::mapped(
        turboXicLoadInputs,
        buildTurboXICLogic
        );
    futures.waitForFinished();

    for (const std::tuple<Err, MzTargetKey, TurboXIC*> &result : futures) {
        e = std::get<0>(result); ree;
        mzTargetKeyVsTurboXicPntr->insert(std::get<1>(result), std::get<2>(result));
    }

    ERR_RETURN

}

Err PythiaDIAFFWorkflowSharedMethods::applyTimsCalibrationEvidencePrefilter(
    const QVector<MsScanInfo> &msScanInfos,
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic &msCalibratomatic,
    const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
    const QMap<MzTargetKey, MsFrame*> &mzTargetKeyVsMsFramePntr,
    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(mzTargetKeyVsTargetDecoyCandidatePointers->isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    e = ErrorUtils::isFalse(mzTargetKeyVsTurboXicPntrs.isEmpty()); ree;
    e = ErrorUtils::isFalse(mzTargetKeyVsMsFramePntr.isEmpty()); ree;

    if (!msCalibratomatic.isInitRT()) {
        ERR_RETURN
    }

    const int candidateBudget = pythiaParameters.timsMainCandidateBudgetPerTargetKey;
    if (candidateBudget <= 0) {
        ERR_RETURN
    }

    for (auto it = mzTargetKeyVsTargetDecoyCandidatePointers->begin();
         it != mzTargetKeyVsTargetDecoyCandidatePointers->end();
         ++it) {

        QVector<TargetDecoyCandidatePair*> &targetDecoyPointers = it.value();
        if (targetDecoyPointers.size() <= candidateBudget) {
            continue;
        }

        const auto turboXicIt = mzTargetKeyVsTurboXicPntrs.constFind(it.key());
        const auto msFrameIt = mzTargetKeyVsMsFramePntr.constFind(it.key());
        if (turboXicIt == mzTargetKeyVsTurboXicPntrs.constEnd()
            || msFrameIt == mzTargetKeyVsMsFramePntr.constEnd()
            || turboXicIt.value() == nullptr
            || msFrameIt.value() == nullptr
            || !turboXicIt.value()->isInit()
            || !msFrameIt.value()->isValid()) {
            continue;
        }

        QVector<TimsEvidenceCandidate> rankedCandidates;
        rankedCandidates.reserve(targetDecoyPointers.size());
        for (TargetDecoyCandidatePair *tdcp : targetDecoyPointers) {
            const double evidenceScore = calculateCheapTimsEvidenceScore(
                *turboXicIt.value(),
                msCalibratomatic,
                pythiaParameters,
                *msFrameIt.value(),
                tdcp
                );
            rankedCandidates.push_back({tdcp, evidenceScore});
        }

        std::sort(
            rankedCandidates.begin(),
            rankedCandidates.end(),
            [](const TimsEvidenceCandidate &left, const TimsEvidenceCandidate &right) {
                if (MathUtils::tSame(left.evidenceScore, right.evidenceScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                    return left.targetDecoyPair->totalFragmentCount() > right.targetDecoyPair->totalFragmentCount();
                }
                return left.evidenceScore > right.evidenceScore;
            }
            );

        QVector<TargetDecoyCandidatePair*> selected;
        selected.reserve(candidateBudget);
        for (const TimsEvidenceCandidate &rankedCandidate : rankedCandidates) {
            if (selected.size() >= candidateBudget) {
                break;
            }
            if (rankedCandidate.evidenceScore <= 0.0) {
                break;
            }
            selected.push_back(rankedCandidate.targetDecoyPair);
        }

        if (!selected.isEmpty()) {
            targetDecoyPointers = selected;
        }
    }

    ERR_RETURN
}

Err PythiaDIAFFWorkflowSharedMethods::buildCandidateScoresPtrs(
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScoresPairsVecBatch,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresPairsVecBatch); ree;

    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> & pr : candidateScoresPairsVecBatch) {
        candidateScoresPntrs->append({&pr.first, &pr.second});
    }

    ERR_RETURN

}

void PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(QVector<CandidateScores*> *candidateScoresPntrs) {
    std::sort(candidateScoresPntrs->begin(), candidateScoresPntrs->end(), [](const CandidateScores *l, const CandidateScores *r){

        if (MathUtils::tSame(l->discriminantScore, r->discriminantScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {

            if (MathUtils::tSame(l->featuresArray[CosineSimSum100], r->featuresArray[CosineSimSum100], S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {

                if (MathUtils::tSame(l->featuresArray[CosineSimSpectrumOverTime], r->featuresArray[CosineSimSpectrumOverTime], S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                    return l->isDecoy > r->isDecoy;
                }

                return l->featuresArray[CosineSimSpectrumOverTime] > r->featuresArray[CosineSimSpectrumOverTime];
            }

            return l->featuresArray[CosineSimSum100] > r->featuresArray[CosineSimSum100];
        }

        return l->discriminantScore > r->discriminantScore;
    });
}

void PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersClassifierScoreAsc(QVector<CandidateScores*> *candidateScoresPntrs) {
    std::sort(
        candidateScoresPntrs->begin(),
        candidateScoresPntrs->end(),
        [](const CandidateScores *l, const CandidateScores *r){
            if (MathUtils::tSame(l->classifierScore, r->classifierScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                return l->discriminantScore > r->discriminantScore;
            }
            return l->classifierScore < r->classifierScore;
        });
}

namespace {

    Err buildProcessBatchPointers(
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScorePairs,
        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy,
        QVector<CandidateScores*> *candidateScoresesPntrs,
        QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> *featuresArrayTargetVsDecoyPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(featuresArrayTargetVsDecoy->isEmpty()); ree;

// #define CHECK_MIN_MAX_MEANS
        QVector<QVector<float>> vecs;
        for (int i = 0; i < featuresArrayTargetVsDecoy->size(); i++) {

            const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &candidatePair = candidateScorePairs.at(i);
            QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &featuresArrayPair = (*featuresArrayTargetVsDecoy)[i];
            featuresArrayTargetVsDecoyPntrs->push_back({&featuresArrayPair.first, &featuresArrayPair.second});

#ifdef CHECK_MIN_MAX_MEANS
            vecs.append({featuresArrayPair.first, featuresArrayPair.second});
#endif
            candidateScoresesPntrs->append({candidatePair.first, candidatePair.second});
        }

#ifdef CHECK_MIN_MAX_MEANS
        qDebug() << "Min max mean check";
        Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        for (int i = 0; i < mat.cols(); i++) {
            const Eigen::VectorX<float> c = mat.col(i);
            qDebug()
            << "Column" << i
            << "min" << c.minCoeff()
            << "max" << c.maxCoeff()
            << "mean" << c.mean()
            << "stdDev" << EigenUtils::calculateStDevOfVector(c)
            ;
        }
#endif

        ERR_RETURN
    }

    Err setCandidateScoresDiscriminantScore(
        const QVector<float> &discriminantScores,
        QVector<CandidateScores*> *candidateScoresesPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(discriminantScores); ree;
        e = ErrorUtils::isEqual(discriminantScores.size(), candidateScoresesPntrs->size()); ree;

        for (int i = 0; i < discriminantScores.size(); i++) {
            CandidateScores* cs = (*candidateScoresesPntrs)[i];
            cs->discriminantScore = discriminantScores.at(i);
        }

        ERR_RETURN
    }

    Err featureSelection(const QVector<FeaturesArray*> &featuresArrayPntrs) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featuresArrayPntrs); ree;

        const Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(featuresArrayPntrs);

        QHash<int, bool> columnsToRemove;

        for (int col1 = 0; col1 < mat.cols(); col1++) {

            const Eigen::VectorX<float> v1 = mat.col(col1);

            for (int col2 = col1 + 1; col2 < mat.cols(); col2++) {

                if (columnsToRemove.value(col1)) {
                    continue;
                }

                const Eigen::VectorX<float> v2 = mat.col(col2);

                double pearsonsCorr;
                e = EigenUtils::pearsonCorrelation(v1, v2, &pearsonsCorr); ree;

                if (std::isnan(pearsonsCorr)) {
                    continue;
                }

                if (constexpr float threshold = 0.9;  pearsonsCorr >= threshold) {
                    qDebug() << "PearsonsCorr" << pearsonsCorr << col1 << col2;
                    columnsToRemove.insert(col2, true);
                }

            }

        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflowSharedMethods::processBatch(
    const QVector<Features> &features,
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScoresPairsVecBatch,
    const PythiaParameters &pythiaParameters,
    QVector<CandidateScores*> *candidateScoresVecBatchPntrs,
    QMap<int, int> *fdrVsCounts,
    QVector<float> *weights,
    bool useMonotonePairQValues
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresPairsVecBatch); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    e = buildCandidateScoresPtrs(
        candidateScoresPairsVecBatch,
        candidateScoresVecBatchPntrs
        ); ree;

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> targetDecoyCandidateScorePairsPntrs;
    e = buildTargetDecoyCandidateScorePairsPntrs(
        candidateScoresPairsVecBatch,
        &targetDecoyCandidateScorePairsPntrs
        ); ree;

    std::sort(
        targetDecoyCandidateScorePairsPntrs.begin(),
        targetDecoyCandidateScorePairsPntrs.end(),
        [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &l, const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &r) {
            return l.first->peptideSequenceWithModsChargeAndTargetKey < r.first->peptideSequenceWithModsChargeAndTargetKey;
        });

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> featuresArrayTargetVsDecoy;
    e = DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
        features,
        targetDecoyCandidateScorePairsPntrs,
        &featuresArrayTargetVsDecoy
        ); ree;

    e = ErrorUtils::isEqual(
        targetDecoyCandidateScorePairsPntrs.size(),
        featuresArrayTargetVsDecoy.size()
        ); ree;

    QVector<CandidateScores*> candidateScoresesPntrs;
    QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> featuresArrayTargetVsDecoyPntrs;
    e = buildProcessBatchPointers(
        targetDecoyCandidateScorePairsPntrs,
        &featuresArrayTargetVsDecoy,
        &candidateScoresesPntrs,
        &featuresArrayTargetVsDecoyPntrs
        ); ree;

    QVector<FeaturesArray*> featuresArrayPntrs;
    for (const QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*> &pr : featuresArrayTargetVsDecoyPntrs) {
        featuresArrayPntrs.append({pr.first, pr.second});
    }

    // e = featureSelection(featuresArrayPntrs); ree;

    e = DiscriminantScoretron::trainLDAClassifier(
            featuresArrayTargetVsDecoyPntrs,
            pythiaParameters.verbosity,
            weights
            ); ree;

    QVector<float> discriminantScores;
    e = DiscriminantScoretron::applyWeights(
        *weights,
        pythiaParameters.threadCount,
        featuresArrayPntrs,
        &discriminantScores
        ); ree;

    e = setCandidateScoresDiscriminantScore(
        discriminantScores,
        &candidateScoresesPntrs
        ); ree;

    //Need to speed this up
    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &targetDecoyCandidateScorePairsPntrs,
        useMonotonePairQValues
        ); ree;

    e = FDRCLassifierNeuralNet::outputFDRResults(
        *candidateScoresVecBatchPntrs,
        pythiaParameters.verbosity,
        fdrVsCounts
        ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflowSharedMethods::buildTargetDecoyCandidateScorePairsPntrs(
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &targetDecoyCandidateScorePairs,
    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidateScorePairs); ree;

    targetDecoyCandidateScorePairsPntrs->clear();

    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &pr : targetDecoyCandidateScorePairs) {
        targetDecoyCandidateScorePairsPntrs->push_back({&pr.first, &pr.second});
    }

    ERR_RETURN
}

namespace {

    void filterTargetDecoyPairPointersByPrecursorMzRange(
        float precursorMzMin,
        float precursorMzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
        ) {

        const auto terminatorLogic = [precursorMzMin, precursorMzMax](const TargetDecoyCandidatePair *tdcp) {
            const float mzPrecursorTargetDecoyPair = tdcp->mz(false);
            return !(precursorMzMin <= mzPrecursorTargetDecoyPair && mzPrecursorTargetDecoyPair <= precursorMzMax);
        };

        const auto terminator = std::remove_if(
            targetDecoyCandidatePairs->begin(),
            targetDecoyCandidatePairs->end(),
            terminatorLogic
            );

        targetDecoyCandidatePairs->erase(terminator, targetDecoyCandidatePairs->end());
    }

} //namespace
Err PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        const PythiaParameters &pythiaParameters,
        const QVector<MsScanInfo> &msScanInfos,
        const MsCalibratomatic *msCalibratomatic,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    mzTargetKeyVsTargetDecoyCandidatePointers->clear();

#define USE_PARALLEL_FILTERING2
#ifdef USE_PARALLEL_FILTERING2
    const auto filterBinder = std::bind(
        filterParallelLogic,
        targetDecoyCandidatePairs,
        pythiaParameters,
        msCalibratomatic,
        std::placeholders::_1
        );

    QFuture<QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>>> futures = QtConcurrent::mapped(
        msScanInfos,
        filterBinder
        );
    futures.waitForFinished();

    for (const QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> &res: futures) {
        e = res.first; ree;
        mzTargetKeyVsTargetDecoyCandidatePointers->insert(res.second.first, res.second.second);
    }
#else
    for (const MsScanInfo &msScanInfo : msScanInfos) {
        QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> result = filterParallelLogic(
            targetDecoyCandidatePairs,
            pythiaParameters,
            msCalibratomatic,
            msScanInfo
            );
        e = result.first; ree;
        mzTargetKeyVsTargetDecoyCandidatePointers->insert(result.second.first, result.second.second);
    }
#endif

    ERR_RETURN
}

QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> PythiaDIAFFWorkflowSharedMethods::filterParallelLogic(
    const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic *msCalibratomatic,
    const MsScanInfo &msScanInfo
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(msScanInfo.msLevel < 2); rree;

    const float precursorMzMin = msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons);
    const float precursorMzMax = msScanInfo.precursorTargetMz + (msScanInfo.isoWindowUpper + pythiaParameters.precursorExtractionWindowThomsons);

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsFiltered;
    targetDecoyCandidatePairsFiltered.reserve(targetDecoyCandidatePairs.size());

    const bool useIonMobilityFilter = msScanInfo.ionMobilityDriftTime > 0.0f;
    for (TargetDecoyCandidatePair *candidate : targetDecoyCandidatePairs) {
        if (candidate == nullptr) {
            continue;
        }

        const float precursorMz = candidate->mz(false);
        if (precursorMz < precursorMzMin || precursorMz > precursorMzMax) {
            continue;
        }

        if (useIonMobilityFilter && !isLibraryIonMobilityInAcquisitionWindow(
            candidate,
            pythiaParameters,
            msCalibratomatic,
            msScanInfo
            )) {
            continue;
        }

        targetDecoyCandidatePairsFiltered.push_back(candidate);
    }

    if (pythiaParameters.verbosity > 0) {
        qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyCandidatePairsFiltered.size() << "targetDecoyPairs found";
    }

    return {e, {msScanInfo.targetKey(), targetDecoyCandidatePairsFiltered}};
}

namespace {

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        std::sort(
            candidateScores->rbegin(),
            candidateScores->rend(),
            [](const CandidateScores *l, const CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
            );

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;

        for (CandidateScores *cs : *candidateScores) {

            const QString key
                = cs->targetDecoyCandidatePair->peptideStringWithMods() + QString::number(cs->targetDecoyCandidatePair->charge());
            if (keyVsCandidatesFoundBest.contains(key)) {
                continue;
            }

            keyVsCandidatesFoundBest.insert(key, cs);
        }

        *candidateScores = keyVsCandidatesFoundBest.values().toVector();
    }
}
Err PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
            const MSLevelEnum &msLevel,
            const QVector<CandidateScores*> &_candidateScores,
            int verbosity,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_candidateScores); ree;

        QVector<CandidateScores*> candidateScoresFiltered = _candidateScores;

        filterDuplicateCandidateScoresByDiscriminantScore(&candidateScoresFiltered);

        if(verbosity > 0) {
            qDebug() << _candidateScores.size() << "Found for recalibartion";
            qDebug() << candidateScoresFiltered.size() << "Found for recalibartion after duplicates filtered";
        }

        constexpr int top6 = 6;
        const auto msCalibrationReaderRowsInsertLogic = [msLevel, top6](CandidateScores *cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
            row.iRTPredicted = cs->targetDecoyCandidatePair->iRt(cs->isDecoy);
            row.scanTime = cs->scanTime;
            row.scanNumber = cs->scanNumber;
            row.driftTime = cs->imDriftTime;
            row.iMPredicted = cs->targetDecoyCandidatePair->iIM();

            if (msLevel == MSLevelEnum::MS2) {

                const QVector<MS2Ion> ms2Ions = cs->isDecoy
                              ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                              : cs->targetDecoyCandidatePair->ms2IonsTarget();

                QVector<float> mzSearchedVals(top6, -1.0f);
                const int maxSize = std::min(top6, ms2Ions.size());

                for (int i = 0; i < maxSize; i++) {
                    mzSearchedVals[i] = ms2Ions.at(i).mz;
                }

                row.mzSearchedVec = mzSearchedVals;
                row.mzFoundMeanVec = cs->featuresArray.mid(MzFoundMean1, top6);
                row.mzFoundStDevVec = cs->featuresArray.mid(MzFoundStDev1, top6);
                row.intensityFoundMaxVec = cs->featuresArray.mid(IntensityFoundMax1, top6);
            }
            else {
                row.mzSearchedVec = {cs->targetDecoyCandidatePair->mz(cs->isDecoy)};
                row.mzFoundMeanVec = {cs->featuresArray[Ms1MzMeanFound100]};
                row.mzFoundStDevVec = {cs->featuresArray[Ms1MzStDevFound100]};
                row.intensityFoundMaxVec = {cs->featuresArray[Ms1IntensityFound100]};
            }

            return row;
        };

        std::transform(
                candidateScoresFiltered.begin(),
                candidateScoresFiltered.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }

namespace {

    std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> buildDiaTargetFramesLoadLogic(
        const MsScanInfo &msi,
        MsReaderPointerAcc *msReaderPointerAcc,
        MsCalibratomatic *msCalibratomatic
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); rtee;

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
        e = msReaderPointerAcc->ptr->getMzTargetScanPoints(msi.targetKey(), &scanNumberVsScanPoints); rtee;

        if (msCalibratomatic->isInitCalMS2()) {
            e = msCalibratomatic->recalibrateScanPoints(
                MSLevelEnum::MS2,
                &scanNumberVsScanPoints
                ); rtee;
        }

        return {e, msi.targetKey(), scanNumberVsScanPoints};
    }

}//namespace
Err PythiaDIAFFWorkflowSharedMethods::buildDiaTargetFrames(
    const QVector<MsScanInfo> &uniqueMsScanInfos,
    MsReaderPointerAcc *msReadPointerAcc,
    MsCalibratomatic *msCalibratomatic,
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames
    ) {

    ERR_INIT

#define RUN_PARALLEL_DIATARGETFRAMES
#ifdef RUN_PARALLEL_DIATARGETFRAMES
QVector<QVector<MsScanInfo>> msScanInfosesTranced;
    e = ParallelUtils::trancheVectorForParallelization(
        uniqueMsScanInfos,
        ParallelUtils::numberOfAvailableSystemProcessors(),
        &msScanInfosesTranced
        );

    const auto trainingLogicBinder = std::bind(
        buildDiaTargetFramesLoadLogic,
        std::placeholders::_1,
        msReadPointerAcc,
        msCalibratomatic
        );

    QFuture<std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>>> futures = QtConcurrent::mapped(
        uniqueMsScanInfos,
        trainingLogicBinder
        );
    futures.waitForFinished();

    for (const std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> &res : futures) {
        e = std::get<0>(res); ree;
        (*diaTargetFrames)[std::get<1>(res)] = std::get<2>(res);
    }
#else
    for (const MsScanInfo& msi : uniqueMsScanInfos) {
        std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> result = buildDiaTargetFramesLoadLogic(
            msi,
            msReadPointerAcc,
            msCalibratomatic
            );
        e = std::get<0>(result); ree;
        (*diaTargetFrames)[std::get<1>(result)] = std::get<2>(result);
    }
#endif

    ERR_RETURN
}

double PythiaDIAFFWorkflowSharedMethods::weightedFDRMean(const QMap<int, int> &fdrVsCounts) {

    QVector<float> weightValues;

    int counter = 0;
    for (int v : fdrVsCounts.values()) {
        weightValues.push_back(v * ((10.0 - counter++)/10.0));
    }

    const double fdrMean = MathUtils::mean(weightValues);

    return fdrMean;
}
