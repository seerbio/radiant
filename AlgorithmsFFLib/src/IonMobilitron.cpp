//
// Created by andrewnichols on 1/7/25.
//

#include "IonMobilitron.h"

#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"

#include "ObjectCSVWriters.h"


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

    Err loadMatrix(
        int column,
        QVector<FeatureFinderHill*> *featureFinderHills,
        Eigen::MatrixX<float> *imMatIntensity,
        Eigen::MatrixX<float> *imMatMz
        ) {

        ERR_INIT

        constexpr int minRows = 1;

        e = ErrorUtils::isAboveThreshold(
            static_cast<int>(imMatIntensity->rows()),
            minRows,
            ErrorUtilsParam::IncludeThreshold
            ); ree;

        e = ErrorUtils::isAboveThreshold(
            static_cast<int>(imMatMz->rows()),
            minRows,
            ErrorUtilsParam::IncludeThreshold
            ); ree;

        const int maxScanNumber = static_cast<int>(imMatIntensity->rows()) - 1;

        for (const FeatureFinderHill *ffh : *featureFinderHills) {

            const QVector<ScanNumberIndex> scanNumberIndexes = ffh->scanNumberIndexes();
            const QVector<double> mzVals = ffh->mzVals();
            const QVector<float> intensityVals = ffh->intensities();

            for (int i = 0; i < scanNumberIndexes.size(); i++) {
                const ScanNumberIndex scanNumberIndex = scanNumberIndexes.at(i);
                const double mzVal = mzVals.at(i);
                const float intensityVal = intensityVals.at(i);

                if (scanNumberIndex > maxScanNumber) {
                    continue;
                }

                imMatIntensity->coeffRef(scanNumberIndex, column) = +intensityVal;
                imMatMz->coeffRef(scanNumberIndex, column) = static_cast<float>(mzVal);
            }
        }

        ERR_RETURN
    }

    Err buildMS1IonMobilityMatriciesAndVectors(
        const Eigen::VectorX<float> &kernel,
        QVector<FeatureFinderHill*> *featureFinderHillsMs1MonoIso,
        QVector<FeatureFinderHill*> *featureFinderHillsMs1Iso1,
        QVector<FeatureFinderHill*> *featureFinderHillsMs1Iso2,
        Eigen::MatrixX<float> *imMatIntensity,
        Eigen::MatrixX<float> *imMatMz,
        Eigen::VectorX<float> *imIntegrationVec
        ) {

        ERR_INIT

        const auto minMaxLogic = [](const FeatureFinderHill *l,const FeatureFinderHill *r) {
            return l->scanNumberIndexMinMax().second < r->scanNumberIndexMinMax().second;
        };

        const FeatureFinderHill *maxIMIndexFFHMono = *std::max_element(
            featureFinderHillsMs1MonoIso->begin(),
            featureFinderHillsMs1MonoIso->end(),
            minMaxLogic
            );

        const IonMobilityIndex maxIMIndex = maxIMIndexFFHMono->scanNumberIndexMinMax().second + 10;

        constexpr int columns = 3;
        imMatIntensity->resize(maxIMIndex, columns);
        imMatIntensity->setZero();

        imMatMz->resize(maxIMIndex, columns);
        imMatMz->setZero();

        e = loadMatrix(
            0,
            featureFinderHillsMs1MonoIso,
            imMatIntensity,
            imMatMz
            ); ree;

        e = loadMatrix(
            1,
            featureFinderHillsMs1Iso1,
            imMatIntensity,
            imMatMz
            ); ree;

        e = loadMatrix(
            2,
            featureFinderHillsMs1Iso2,
            imMatIntensity,
            imMatMz
            ); ree;

        *imMatIntensity = EigenKernelUtils::applyKernelToEachColumnInMatrix(
            *imMatIntensity,
            kernel
            );

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;
        Eigen::MatrixX<float> imMatIntensityDigital = *imMatIntensity;
        imMatIntensityDigital
            = (imMatIntensityDigital.array() > intensityThresholdVal).select(countValue, imMatIntensityDigital);
        EigenUtils::thresholdMatrix(0.0f, &imMatIntensityDigital);

        *imIntegrationVec = imMatIntensityDigital.rowwise().sum();
        *imIntegrationVec = EigenKernelUtils::convolveVectorWithKernel(*imIntegrationVec, kernel);
        *imIntegrationVec = imIntegrationVec->array() * imMatIntensity->array().rowwise().sum();

        ERR_RETURN;
    }

    float calculateWeightedMeanOfIndexes(
        const Eigen::VectorX<float> &imIntegrationVec,
        const PeakIntegrationIndexes &imPeakIntegrationIndexes
        ) {

        float numerator = 0.0f;
        float denominator = 0.0f;

        const int peakLen = imPeakIntegrationIndexes.first + (imPeakIntegrationIndexes.second - imPeakIntegrationIndexes.first + 1);

        for (int i = imPeakIntegrationIndexes.first; i < peakLen; i++) {

            const float intensity = imIntegrationVec.coeff(i);

            if (MathUtils::tZero(intensity)) {
                continue;
            }

            numerator += i * intensity;
            denominator += intensity;
        }

        if (MathUtils::tZero(denominator)) {
            return -1;
        }

        return static_cast<int>(std::round(numerator / denominator));
    }

    Err processMS1IonMobilityIntegrationVec(
        const Eigen::VectorX<float> &imIntegrationVec,
        const Eigen::VectorX<float> &kernel,
        int *imPeakCount,
        float *imArea,
        float *imAreaRatioTop2,
        int *imIndexMax,
        PeakIntegrationIndexes *imPeakIntegrationIndexes
        ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(
            static_cast<int>(imIntegrationVec.rows()),
            0,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

        const Eigen::VectorX<float> imIntegrationVecSmoothed = EigenKernelUtils::convolveVectorWithKernel(
            imIntegrationVec,
            kernel
            );

        constexpr float stopThresholdFractionMS2 = 0.2;
        constexpr float peakDifferenceFractionThreshold = 0.4;
        QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
        e = EigenUtils::simpleIntegrator(
            imIntegrationVecSmoothed,
            stopThresholdFractionMS2,
            peakDifferenceFractionThreshold,
            &peakIntegrationsVsIntensities
            ); ree;

        *imPeakCount = peakIntegrationsVsIntensities.size();

        using CONT = QPair<PeakIntegrationIndexes, Intensity>;
        std::sort(
            peakIntegrationsVsIntensities.rbegin(),
            peakIntegrationsVsIntensities.rend(),
            [](const CONT &l, const CONT &r) {return l.second < r.second;}
            );

        QVector<float> trapAreas;
        for (int i = 0; i < std::min(peakIntegrationsVsIntensities.size(), 2); i++) {

            const QPair<PeakIntegrationIndexes, Intensity> pr = peakIntegrationsVsIntensities.at(i);
            const Eigen::VectorX<float> imIntegrationVecSmoothedPeak = imIntegrationVecSmoothed.segment(
                pr.first.first,
                pr.first.second - pr.first.first + 1
            );

            trapAreas.push_back(EigenUtils::calculateTrapizoidalArea(imIntegrationVecSmoothedPeak));
        }

        if (trapAreas.isEmpty()) {
            *imArea = 0.0;
            *imAreaRatioTop2 = 0.0f;
            *imPeakIntegrationIndexes = {-1, -1};
            *imIndexMax = -1;
        }
        else if (trapAreas.size() < 2) {
            *imArea = trapAreas.front();
            *imAreaRatioTop2 = 0.0f;
            *imPeakIntegrationIndexes = peakIntegrationsVsIntensities.at(0).first;
            *imIndexMax = calculateWeightedMeanOfIndexes(imIntegrationVec, *imPeakIntegrationIndexes);
        }
        else {

            const int trapAreaIndexMax = MathUtils::findMaxIndexInVector(trapAreas);
            const int trapAreaIndexMin = trapAreaIndexMax == 1 ? 0 : 1;

            *imArea = trapAreas.at(trapAreaIndexMax);
            *imAreaRatioTop2 = trapAreas.at(trapAreaIndexMax) / trapAreas.at(trapAreaIndexMin);
            *imPeakIntegrationIndexes = peakIntegrationsVsIntensities.at(trapAreaIndexMax).first;
            *imIndexMax = calculateWeightedMeanOfIndexes(imIntegrationVec, *imPeakIntegrationIndexes);
        }

        ERR_RETURN;
    }

    Err assignIonMobilityLogic(
            float ppmTolerance,
            const QVector<CandidateScores*> &candidateScoresPntrs,
            const Eigen::MatrixX<float> &kernel,
            QMap<ScanNumber, FeatureFinderHillBuilder*> *scanNumberVsFeatureFinderHillBuildersPntrsTIMS
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;
        e = ErrorUtils::isFalse(scanNumberVsFeatureFinderHillBuildersPntrsTIMS->isEmpty()); ree;

        const QVector<ScanNumber> scanNumbers = scanNumberVsFeatureFinderHillBuildersPntrsTIMS->keys().toVector();

        for (CandidateScores *cs : candidateScoresPntrs) {

            if (cs->scanNumber < 0) {
                continue;
            }

            FrameIndex ms1ScanNumberClosest = scanNumbers.at(MathUtils::closest(scanNumbers, cs->scanNumber));
            ms1ScanNumberClosest = ms1ScanNumberClosest > cs->scanNumber
                                 ? scanNumbers.at(MathUtils::closest(scanNumbers, cs->scanNumber) - 1)
                                 : ms1ScanNumberClosest;

            e = ErrorUtils::contains(ms1ScanNumberClosest, *scanNumberVsFeatureFinderHillBuildersPntrsTIMS); ree;

            const FeatureFinderHillBuilder *ffhb = scanNumberVsFeatureFinderHillBuildersPntrsTIMS->value(ms1ScanNumberClosest);

            const float mzMonoIso = cs->targetDecoyCandidatePair->mz(cs->isDecoy);
            const float massToleranceMonoIso = MathUtils::calculatePPM(mzMonoIso, ppmTolerance);
            QVector<FeatureFinderHill*> featureFinderHillsMs1MonoIso;
            e = ffhb->getHills(
                static_cast<double>(mzMonoIso - massToleranceMonoIso),
                static_cast<double>(mzMonoIso + massToleranceMonoIso),
                &featureFinderHillsMs1MonoIso
                ); ree;

            const float isoDistance
                = S_GLOBAL_SETTINGS.ISO_DIFF / static_cast<float>(cs->targetDecoyCandidatePair->charge());

            const float mzIso1 = mzMonoIso + isoDistance;
            const float massToleranceIso1 = MathUtils::calculatePPM(mzIso1, ppmTolerance);
            QVector<FeatureFinderHill*> featureFinderHillsMs1Iso1;
            e = ffhb->getHills(
                static_cast<double>(mzIso1 - massToleranceIso1),
                static_cast<double>(mzIso1 + massToleranceIso1),
                &featureFinderHillsMs1Iso1
                ); ree;

            const float mzIso2 = mzMonoIso + (isoDistance * 2);
            const float massToleranceIso2 = MathUtils::calculatePPM(mzIso2, ppmTolerance);
            QVector<FeatureFinderHill*> featureFinderHillsMs1Iso2;
            e = ffhb->getHills(
                static_cast<double>(mzIso2 - massToleranceIso2),
                static_cast<double>(mzIso2 + massToleranceIso2),
                &featureFinderHillsMs1Iso2
                ); ree;

            if (featureFinderHillsMs1MonoIso.isEmpty()) {
                continue;
            }

            Eigen::MatrixX<float> imMatIntensity;
            Eigen::MatrixX<float> imMatMz;
            Eigen::VectorX<float> imIntegrationVec;
            e = buildMS1IonMobilityMatriciesAndVectors(
                kernel,
                &featureFinderHillsMs1MonoIso,
                &featureFinderHillsMs1Iso1,
                &featureFinderHillsMs1Iso2,
                &imMatIntensity,
                &imMatMz,
                &imIntegrationVec
                ); ree;

            int imPeakCount;
            float imArea;
            float imAreaRatioTop2;
            int imIndexMax;
            PeakIntegrationIndexes imPeakIntegrationIndexes;
            e = processMS1IonMobilityIntegrationVec(
                imIntegrationVec,
                kernel,
                &imPeakCount,
                &imArea,
                &imAreaRatioTop2,
                &imIndexMax,
                &imPeakIntegrationIndexes
                ); ree;

            const int imPeakLen = imPeakIntegrationIndexes.second - imPeakIntegrationIndexes.first + 1;

            // cs->featuresArray[Features::ImPeakCount] = static_cast<float>(imPeakCount);
            // cs->featuresArray[Features::ImPeakLength] = imPeakLen < 2 ? -1 : static_cast<float>(imPeakLen);
            // cs->featuresArray[Features::ImAreaLog10] = std::log10(imArea);
            // cs->featuresArray[Features::ImAreaRatioTop2] = imAreaRatioTop2;
            cs->ionMobilityIndex = imIndexMax;
            cs->ionMobilityIndexStart = imPeakIntegrationIndexes.first;
            cs->ionMobilityIndexEnd = imPeakIntegrationIndexes.second;
        }

        ERR_RETURN
    }

}//namespace
Err IonMobilitron::assignIonMobilityValues(
    const PythiaParameters &pythiaParameters,
    const QVector<CandidateScores*> &_candidateScoresPntrs,
    QMap<ScanNumber, FeatureFinderHillBuilder*> *scanNumberVsFeatureFinderHillBuildersPntrsTIMS
    ) {

    ERR_INIT

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "IonMobilitron::assignIonMobilityValues";

    e = ErrorUtils::isNotEmpty(_candidateScoresPntrs); ree;

    QVector<CandidateScores*> candidateScoresPntrs = _candidateScoresPntrs;

    constexpr int filterLength = 3;
    constexpr int order = 1;
    constexpr int derivative = 0;
    constexpr int rate = 1;

    Eigen::MatrixX<float> kernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        filterLength,
        order,
        derivative,
        rate,
        &kernel
        ); ree;
    const Eigen::VectorX<float> kernelVec(kernel);

    std::sort(
        candidateScoresPntrs.begin(),
        candidateScoresPntrs.end(),
        [](const CandidateScores* l, const CandidateScores* r) {return l->scanNumber < r->scanNumber;}
        );

    QVector<QVector<CandidateScores*>> candidateScoresPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelizationInOrder(
        candidateScoresPntrs,
        pythiaParameters.threadCount,
        0,
        &candidateScoresPntrsTranched
        ); ree;

#define ASSIGN_ION_MOBILITY_PARALLEL
#ifdef ASSIGN_ION_MOBILITY_PARALLEL
    const auto assignMobilityBinder = std::bind(
            assignIonMobilityLogic,
            pythiaParameters.ms2ExtractionWidthPPM,
            std::placeholders::_1,
            kernelVec,
            scanNumberVsFeatureFinderHillBuildersPntrsTIMS
            );

    QFuture<Err> futures = QtConcurrent::mapped(
        candidateScoresPntrsTranched,
        assignMobilityBinder
        );
    futures.waitForFinished();
#else
    for (QVector<CandidateScores*> &csTranche : candidateScoresPntrsTranched) {
        e = assignIonMobilityLogic(
            pythiaParameters.ms2ExtractionWidthPPM,
            csTranche,
            kernel,
            scanNumberVsFeatureFinderHillBuildersPntrsTIMS
            ); ree;
    }
#endif

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "IonMobilitron::assignIonMobilityValues finished";

    ERR_RETURN
}
