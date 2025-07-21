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

    struct IonMobilityInput {
        CandidateScores *candidateScores = nullptr;
        QMap<FrameNumberTIMS, Ms1FrameTIMS> *ms1Frames = nullptr;
        MsReaderPointerAcc *msReaderPointerAcc = nullptr;
        float ppmTol = -1.0;
        QVector<ScanNumber> scanNumbers;
    };

    Err assignIonMobilityIndexesToCandidateScoresLogic(const IonMobilityInput &imi) {

        ERR_INIT

        FrameIndex ms1ScanNumberClosest = imi.scanNumbers.at(MathUtils::closest(imi.scanNumbers, imi.candidateScores->scanNumber));
        ms1ScanNumberClosest = ms1ScanNumberClosest > imi.candidateScores->scanNumber
                             ? imi.scanNumbers.at(MathUtils::closest(imi.scanNumbers, imi.candidateScores->scanNumber) - 1)
                             : ms1ScanNumberClosest;

        e = ErrorUtils::contains(ms1ScanNumberClosest, *imi.ms1Frames); ree;
        Ms1FrameTIMS ms1FrameTims = imi.ms1Frames->value(ms1ScanNumberClosest);

        const std::tuple<Err, float, float> result = IonMobilitron::findEmpericalIonMobilityDriftTime(
            imi.candidateScores,
            &ms1FrameTims,
            imi.msReaderPointerAcc,
            imi.ppmTol
            );
        e = std::get<0>(result); ree;

        imi.candidateScores->imDriftTime = std::get<1>(result);
        imi.candidateScores->featuresArray[Ms1IntensityFoundApex100IM] = std::get<2>(result);

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

    QVector<IonMobilityInput> ionMobilityInputs;
    for (CandidateScores *cs : candidateScoresVecBatchPntrs) {

        const float cosineSim100MS1 = cs->featuresArray[CosineSim100MS1];
        if (constexpr float cutOff = 0.9; cosineSim100MS1 < cutOff) {
            continue;
        }

        IonMobilityInput imi;
        imi.candidateScores = cs;
        imi.ms1Frames = ms1Frames;
        imi.ppmTol = ppmTol;
        imi.scanNumbers = scanNumbers;
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

    Err getXICPoints(
        Ms1FrameTIMS *ms1FrameTims,
        float mz,
        int charge, float extractionPPMTol, QVector<XICPoints> *xicPointses
    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(ms1FrameTims->isEmpty()); ree;

        QMap<ScanNumber, ScanPoints*> ms1FrameTimsPntrs;
        for (auto it = ms1FrameTims->begin(); it != ms1FrameTims->end(); ++it) {
            ms1FrameTimsPntrs.insert(it.key(), &it.value());
        }

        TurboXIC turboXic;
        e = turboXic.init(ms1FrameTimsPntrs); ree;

        const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / charge;
        const float mzMono = mz;

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
        Eigen::MatrixX<float> *xicMat,
        Eigen::VectorX<float> *integrationVector
        ) {

        ERR_INIT

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

        const Eigen::VectorX<float> sgKernelVec = sgKernel;
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

    QVector<XICPoints> xicPointses;
    e = getXICPoints(
        ms1FrameTims,
        candidateScores->mz(),
        candidateScores->charge(),
        extractionPPMTol, &xicPointses
    ); rtee;

    const XICPoints &xicPointsMonoIsotope = xicPointses.front();
    if (xicPointsMonoIsotope.empty()) {
        return {e, -1.0f, -1.0f};
    }

    Eigen::VectorX<float> integrationVector;
    Eigen::MatrixX<float> xicMat;
    e = buildCombinedXIC(
        xicPointses,
        &xicMat,
        &integrationVector
        ); rtee;

    constexpr float stopThresholdFractionMS2 = 0.2;
    constexpr float peakDifferenceFractionThreshold = 0.25;

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        integrationVector,
        stopThresholdFractionMS2,
        peakDifferenceFractionThreshold,
        &peakIntegrationsVsIntensities
        ); rtee;

    const QPair<PeakIntegrationIndexes, Intensity> &prImLimits = peakIntegrationsVsIntensities.front();

    const XICPoint xicPointApex = *std::max_element(
        xicPointsMonoIsotope.begin(),
        xicPointsMonoIsotope.end(),
        [](const XICPoint &l, const XICPoint &r) {return l.intensity < r.intensity;}
        );

    double driftTime;
    e = msReaderPointerAcc->ptr->driftTimeFromIonMobilityIndex(
        xicPointApex.scanNumber,
        &driftTime
        ); rtee;

    return {e, static_cast<float>(driftTime), xicPointApex.intensity};
}
