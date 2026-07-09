//
// Created by andrewnichols on 1/7/25.
//

#include "IonMobilitron.h"

#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsReaderBase.h"
#include "ObjectCSVWriters.h"
#include "TurboXIC.h"

#include <memory>
#include <mutex>


Err IonMobilitron::init(const QVector<QPair<IMPredicted, IMEmpirical>> &imPredVsImEmpValuesSortedDiscScoreHiLo) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(imPredVsImEmpValuesSortedDiscScoreHiLo); ree;

    const QVector<int> heads = {400, 350, 300, 250, 200, 150, 100, 50};

    QVector<QPair<IMPredicted, IMEmpirical>> imPredVsImEmpValuesSortedDiscScoreHiLoCopy = imPredVsImEmpValuesSortedDiscScoreHiLo;

    double bestPValue = std::numeric_limits<double>::max();
    int bestHead = 0;
    for (int head : heads) {
        imPredVsImEmpValuesSortedDiscScoreHiLoCopy.resize(head);

        double slope;
        double intercept;
        double stdError;
        double pValue;
        e = EigenUtils::linearRegression(
            imPredVsImEmpValuesSortedDiscScoreHiLoCopy,
            &slope,
            &intercept,
            &stdError,
            &pValue
            ); ree;

        qDebug()
        << qPrintable(S_GLOBAL_TIMER.elapsed())
        << "IonMobilitron init optimiztation:"
        << "head count" << head
        << "| slope" << slope
        << "| intercept" << intercept
        << "| std_error" << stdError
        << "| p_value" << pValue;

        if (bestPValue > stdError && slope < 0.0) {
            e = m_iIMtoIonMobilityIndexMapper.init(imPredVsImEmpValuesSortedDiscScoreHiLoCopy); ree;
            bestPValue = stdError;
            bestHead = head;
            // for (const QPair<IMPredicted, IMEmpirical> &pr : imPredVsImEmpValuesSortedDiscScoreHiLoCopy) {
            //     qDebug() << "(" << pr.first << ", " << pr.second << ")";
            //
            // }
        }
    }

    qDebug()
    << qPrintable(S_GLOBAL_TIMER.elapsed())
    << "Head count of" << bestHead << "chosen for IonMobilitron init";

    ERR_RETURN
}

Err IonMobilitron::predictIonMobilityIndex(
    float iIM,
    int *predictedIonMobilityIndex
    ) const {

    ERR_INIT

    double predictedScanTimeDouble;
    e = m_iIMtoIonMobilityIndexMapper.predictY(static_cast<double>(iIM), &predictedScanTimeDouble); ree;
    *predictedIonMobilityIndex = static_cast<int>(predictedScanTimeDouble);

    ERR_RETURN
}

namespace {

    struct IonMobilityFrameCache {
        const Ms1FrameTIMS *ms1FrameTims = nullptr;
        QMap<ScanNumber, ScanPoints*> ms1FrameTimsPntrs;
        TurboXIC turboXic;
        std::once_flag initOnce;
        Err initErr = eNoError;
        bool isInit = false;
    };

    struct IonMobilityInput {
        CandidateScores *candidateScores = nullptr;
        std::shared_ptr<IonMobilityFrameCache> frameCache;
        MsReaderPointerAcc *msReaderPointerAcc = nullptr;
        std::shared_ptr<const Eigen::VectorX<float>> sgKernelVec;
        float ppmTol = -1.0;
    };

    Err buildSavitzkyGolayKernelVec(Eigen::VectorX<float> *sgKernelVec) {

        ERR_INIT

        e = ErrorUtils::isFalse(sgKernelVec == nullptr); ree;

        constexpr int filterLength = 3;
        constexpr int order = 1;
        constexpr int derivative = 0;
        constexpr int rate = 1;

        Eigen::MatrixX<float> sgKernel;
        e = EigenKernelUtils::buildSavitzkyGolayKernel(
            filterLength,
            order,
            derivative,
            rate,
            &sgKernel
            ); ree;

        *sgKernelVec = sgKernel;

        ERR_RETURN
    }

    bool nearestPrecedingScanNumber(
        const QVector<ScanNumber> &scanNumbers,
        ScanNumber scanNumber,
        ScanNumber *nearestScanNumber
        ) {

        if (nearestScanNumber == nullptr || scanNumbers.isEmpty()) {
            return false;
        }

        int closestIndex = MathUtils::closest(scanNumbers, scanNumber);
        if (closestIndex < 0 || closestIndex >= scanNumbers.size()) {
            return false;
        }

        ScanNumber scanNumberClosest = scanNumbers.at(closestIndex);
        if (scanNumberClosest > scanNumber) {
            if (closestIndex == 0) {
                return false;
            }
            scanNumberClosest = scanNumbers.at(closestIndex - 1);
        }

        *nearestScanNumber = scanNumberClosest;
        return true;
    }

    Err buildIonMobilityFrameCache(IonMobilityFrameCache *frameCache) {

        ERR_INIT

        e = ErrorUtils::isFalse(frameCache == nullptr); ree;
        e = ErrorUtils::isFalse(frameCache->ms1FrameTims == nullptr); ree;
        e = ErrorUtils::isFalse(frameCache->ms1FrameTims->isEmpty()); ree;

        frameCache->ms1FrameTimsPntrs.clear();
        for (auto it = frameCache->ms1FrameTims->constBegin(); it != frameCache->ms1FrameTims->constEnd(); ++it) {
            frameCache->ms1FrameTimsPntrs.insert(it.key(), const_cast<ScanPoints*>(&it.value()));
        }

        e = frameCache->turboXic.init(frameCache->ms1FrameTimsPntrs); ree;
        frameCache->isInit = frameCache->turboXic.isInit();

        ERR_RETURN
    }

    Err ensureIonMobilityFrameCacheInit(IonMobilityFrameCache *frameCache) {

        if (frameCache == nullptr) {
            return eError;
        }

        std::call_once(frameCache->initOnce, [frameCache]() {
            frameCache->initErr = buildIonMobilityFrameCache(frameCache);
        });

        return frameCache->initErr;
    }

    Err getXICPoints(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const TurboXIC &turboXic,
        float extractionPPMTol,
        QVector<XICPoints> *xicPointses
        ) {

        ERR_INIT

        const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / targetDecoyCandidatePair->charge();
        const float mzMono = targetDecoyCandidatePair->mz(false);

        for (int i = 0; i < 2; i++) {
            const float mzExtract = mzMono + (i * chargeDistance);
            const float massTol = MathUtils::calculatePPM(mzExtract, extractionPPMTol);

            const XICPoints xicPoints = turboXic.extractPointsXIC(
                mzExtract - massTol,
                mzExtract + massTol
                );
            xicPointses->push_back(xicPoints);
        }

        ERR_RETURN
    }

    Err buildCombinedXIC(
        const QVector<XICPoints> &xicPointses,
        const Eigen::VectorX<float> &sgKernelVec,
        Eigen::MatrixX<float> *xicMat,
        Eigen::VectorX<float> *integrationVector
        );

    Err findEmpericalIonMobilityDriftTimeLogic(
        CandidateScores *candidateScores,
        IonMobilityFrameCache &frameCache,
        MsReaderPointerAcc *msReaderPointerAcc,
        float extractionPPMTol,
        const Eigen::VectorX<float> &sgKernelVec,
        float *empiricalDriftTime,
        float *apexIntensity
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(frameCache.ms1FrameTims == nullptr); ree;
        e = ErrorUtils::isFalse(frameCache.ms1FrameTims->isEmpty()); ree;
        e = ErrorUtils::isTrue(frameCache.turboXic.isInit()); ree;
        e = ErrorUtils::isFalse(empiricalDriftTime == nullptr); ree;
        e = ErrorUtils::isFalse(apexIntensity == nullptr); ree;

        QVector<XICPoints> xicPointses;
        e = getXICPoints(
            candidateScores->targetDecoyCandidatePair,
            frameCache.turboXic,
            extractionPPMTol,
            &xicPointses
            ); ree;

        const XICPoints &xicPointsMonoIsotope = xicPointses.front();
        if (xicPointsMonoIsotope.empty()) {
            *empiricalDriftTime = -1.0f;
            *apexIntensity = -1.0f;
            ERR_RETURN
        }

        Eigen::VectorX<float> integrationVector;
        Eigen::MatrixX<float> xicMat;
        e = buildCombinedXIC(
            xicPointses,
            sgKernelVec,
            &xicMat,
            &integrationVector
            ); ree;

        constexpr float stopThresholdFractionMS2 = 0.2;
        constexpr float peakDifferenceFractionThreshold = 0.25;

        QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
        e = EigenUtils::simpleIntegrator(
            integrationVector,
            stopThresholdFractionMS2,
            peakDifferenceFractionThreshold,
            &peakIntegrationsVsIntensities
            ); ree;

        const QPair<PeakIntegrationIndexes, Intensity> &prImLimits = peakIntegrationsVsIntensities.front();
        Q_UNUSED(prImLimits)

        const XICPoint xicPointApex = *std::max_element(
            xicPointsMonoIsotope.begin(),
            xicPointsMonoIsotope.end(),
            [](const XICPoint &l, const XICPoint &r) {return l.intensity < r.intensity;}
            );

        double driftTime;
        e = msReaderPointerAcc->ptr->driftTimeFromIonMobilityIndex(
            xicPointApex.scanNumber,
            &driftTime
            ); ree;

        *empiricalDriftTime = static_cast<float>(driftTime);
        *apexIntensity = xicPointApex.intensity;

        ERR_RETURN
    }

    Err assignIonMobilityIndexesToCandidateScoresLogic(const IonMobilityInput &imi) {

        ERR_INIT

        e = ErrorUtils::isFalse(imi.frameCache.get() == nullptr); ree;
        e = ErrorUtils::isFalse(imi.msReaderPointerAcc == nullptr); ree;
        e = ErrorUtils::isFalse(imi.sgKernelVec.get() == nullptr); ree;

        e = ensureIonMobilityFrameCacheInit(imi.frameCache.get()); ree;
        if (!imi.frameCache->isInit) {
            ERR_RETURN
        }

        float empiricalDriftTime = -1.0f;
        float apexIntensity = -1.0f;
        e = findEmpericalIonMobilityDriftTimeLogic(
            imi.candidateScores,
            *imi.frameCache,
            imi.msReaderPointerAcc,
            imi.ppmTol,
            *imi.sgKernelVec,
            &empiricalDriftTime,
            &apexIntensity
            ); ree;

        imi.candidateScores->imDriftTime = empiricalDriftTime;
        imi.candidateScores->featuresArray[Ms1IntensityFoundApex100IM] = apexIntensity;

        ERR_RETURN
    }

}//namespace
Err IonMobilitron::assignIonMobilityIndexesToCandidateScores(
    const QVector<CandidateScores*> &candidateScoresVecBatchPntrs,
    float ppmTol,
    MsReaderPointerAcc *msReaderPointerAcc
    ) {

    ERR_INIT

    QMap<FrameNumberTIMS, Ms1FrameTIMS> *ms1Frames = msReaderPointerAcc->ptr->frameNumberVsMS1FrameTIMSPntr();
    const QMap<ScanNumber, MsScanInfo> ms1ScanInfos = msReaderPointerAcc->ptr->getMsScanInfos(1);
    const QVector<ScanNumber> scanNumbers = ms1ScanInfos.keys().toVector();
    auto sgKernelVec = std::make_shared<Eigen::VectorX<float>>();
    e = buildSavitzkyGolayKernelVec(sgKernelVec.get()); ree;

    QMap<ScanNumber, std::shared_ptr<IonMobilityFrameCache>> frameCacheByScanNumber;
    QVector<IonMobilityInput> ionMobilityInputs;
    ionMobilityInputs.reserve(candidateScoresVecBatchPntrs.size());
    for (CandidateScores *cs : candidateScoresVecBatchPntrs) {

        const float cosineSim100MS1 = cs->featuresArray[CosineSim100MS1];
        if (constexpr float cutOff = 0.9; cosineSim100MS1 < cutOff) {
            continue;
        }

        ScanNumber ms1ScanNumberClosest = -1;
        if (!nearestPrecedingScanNumber(
            scanNumbers,
            cs->scanNumber,
            &ms1ScanNumberClosest
            )) {
            continue;
        }

        const auto ms1FrameIt = ms1Frames->constFind(ms1ScanNumberClosest);
        if (ms1FrameIt == ms1Frames->constEnd() || ms1FrameIt.value().isEmpty()) {
            continue;
        }

        std::shared_ptr<IonMobilityFrameCache> frameCache = frameCacheByScanNumber.value(ms1ScanNumberClosest);
        if (!frameCache) {
            frameCache = std::make_shared<IonMobilityFrameCache>();
            frameCache->ms1FrameTims = &ms1FrameIt.value();
            frameCacheByScanNumber.insert(ms1ScanNumberClosest, frameCache);
        }

        IonMobilityInput imi;
        imi.candidateScores = cs;
        imi.frameCache = frameCache;
        imi.ppmTol = ppmTol;
        imi.sgKernelVec = sgKernelVec;
        imi.msReaderPointerAcc = msReaderPointerAcc;
        ionMobilityInputs.push_back(imi);
    }

#define PARALLEL_ION_MOBILITY
#ifdef PARALLEL_ION_MOBILITY
    QFuture<Err> futures = QtConcurrent::mapped(
        ionMobilityInputs,
        assignIonMobilityIndexesToCandidateScoresLogic
        );
    futures.waitForFinished();
#else
    for (IonMobilityInput &imi : ionMobilityInputs) {
        e = assignIonMobilityIndexesToCandidateScoresLogic(imi); ree;
    }
#endif//PARALLEL_ION_MOBILITY

    ERR_RETURN
}

namespace {
    Err buildCombinedXIC(
        const QVector<XICPoints> &xicPointses,
        const Eigen::VectorX<float> &sgKernelVec,
        Eigen::MatrixX<float> *xicMat,
        Eigen::VectorX<float> *integrationVector
        ) {

        ERR_INIT

        const auto sortLogic
            = [](const XICPoint &l, const XICPoint &r) {return l.scanNumber < r.scanNumber;};

        IonMobilityIndex ionMobilityIndexMax = 0;

        for (const XICPoints &xicPoints : xicPointses ) {

            if (xicPoints.empty()) {
                continue;
            }

            const IonMobilityIndex maxIonMobilityIndex = std::max_element(
                xicPoints.begin(),
                xicPoints.end(),
                sortLogic
                )->scanNumber;

            ionMobilityIndexMax = std::max(ionMobilityIndexMax, maxIonMobilityIndex);
        }

        if (ionMobilityIndexMax < 1) {
            ERR_RETURN
        }

        constexpr int buffer = 10;
        integrationVector->resize(ionMobilityIndexMax + buffer);
        integrationVector->setZero();
        xicMat->resize(ionMobilityIndexMax + buffer, 3);
        xicMat->setZero();

        int column = 0;
        for (const XICPoints &xicPoints : xicPointses) {
            for (const XICPoint &xicPoint : xicPoints) {

                e = ErrorUtils::isBelowThreshold(
                    xicPoint.scanNumber,
                    ionMobilityIndexMax + buffer,
                    ErrorUtilsParam::ExcludeThreshold
                    ); ree;

                integrationVector->coeffRef(xicPoint.scanNumber) += 1.0f;
                xicMat->coeffRef(xicPoint.scanNumber, column) += xicPoint.intensity;
            }
            column++;
        }

        const Eigen::VectorX<float> isotopeIntenstiesSummed = xicMat->rowwise().sum();

        *integrationVector = integrationVector->array() * isotopeIntenstiesSummed.array();

        *integrationVector = EigenKernelUtils::convolveVectorWithKernel(
            *integrationVector,
            sgKernelVec
            );

        ERR_RETURN
    }

}//namespace
std::tuple<Err, float, float> IonMobilitron::findEmpericalIonMobilityDriftTime(
    CandidateScores *candidateScores,
    Ms1FrameTIMS *ms1FrameTims,
    MsReaderPointerAcc *msReaderPointerAcc,
    float extractionPPMTol
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(ms1FrameTims->isEmpty()); rtee;

    IonMobilityFrameCache frameCache;
    frameCache.ms1FrameTims = ms1FrameTims;
    e = buildIonMobilityFrameCache(&frameCache); rtee;

    Eigen::VectorX<float> sgKernelVec;
    e = buildSavitzkyGolayKernelVec(&sgKernelVec); rtee;

    float empiricalDriftTime = -1.0f;
    float apexIntensity = -1.0f;
    e = findEmpericalIonMobilityDriftTimeLogic(
        candidateScores,
        frameCache,
        msReaderPointerAcc,
        extractionPPMTol,
        sgKernelVec,
        &empiricalDriftTime,
        &apexIntensity
        ); rtee;

    return {e, empiricalDriftTime, apexIntensity};
}
