//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron2.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "IsotopicDistributionBuilder.h"
#include "ScoreOverseer.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"


class Q_DECL_HIDDEN CandidateScorertron::Private {
public:

    using NominalMass = int;
    using AveragineVec = QVector<float>;

    Private();
    ~Private();

    Err init();

    Eigen::VectorX<float> m_kernelIntegration;
    Eigen::VectorX<float> m_kernelMs2;

};

CandidateScorertron::Private::Private() {}
CandidateScorertron::Private::~Private() {}

Err CandidateScorertron::Private::init() {

    ERR_INIT

    constexpr int filterLengthMs2 = 5;
    constexpr int filterLengthIntegration = 5;
    constexpr int order = 2;
    constexpr int derivative = 0;
    constexpr int rate = 1;

    Eigen::MatrixX<float> kernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        filterLengthMs2,
        order,
        derivative,
        rate,
        &kernel
        ); ree;
    const Eigen::VectorX<float> kernelVec(kernel);
    m_kernelMs2 = kernelVec;

    Eigen::MatrixX<float> kernelIntegration;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        filterLengthIntegration,
        order,
        derivative,
        rate,
        &kernelIntegration
        ); ree;
    const Eigen::VectorX<float> kernelIntegrationVec(kernelIntegration);
    m_kernelIntegration = kernelIntegrationVec;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_xicPeakManager(nullptr)
, m_msFrameMzTarget(nullptr)
, d_ptr(new Private())
{}

CandidateScorertron::~CandidateScorertron() {}

Err CandidateScorertron::init(
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic &msCalibratomatic,
    int topNMS2Ions,
    XICPeakManager *xicPeakManager,
    MsFrame *msFrameMzTarget
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
    e = ErrorUtils::isAboveThreshold(
        topNMS2Ions,
        S_GLOBAL_SETTINGS.MIN_MS2_IONS,
        ErrorUtilsParam::IncludeThreshold)
    ; ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_xicPeakManager = xicPeakManager;
    m_msFrameMzTarget = msFrameMzTarget;

    if (msCalibratomatic.isInit()) {
        m_msCalibratomatic = msCalibratomatic;
    }

    e = d_ptr->init(); ree;

    ERR_RETURN
}

class MatriciesAndVecs {

public:

    Eigen::MatrixX<float> intensityMatrix100;
    Eigen::MatrixX<float> intensityMatrix100Shadow;
    Eigen::MatrixX<float> intensityMatrix45;
    Eigen::MatrixX<float> intensityMatrix20;

    Eigen::MatrixX<float> mzMatrix100;
    Eigen::MatrixX<float> mzMatrix100Shadow;
    Eigen::MatrixX<float> mzMatrix45;
    Eigen::MatrixX<float> mzMatrix20;

    Eigen::VectorX<float> integrationVec;

    [[nodiscard]] bool intensityMatriciesAreValid() const {
        return intensityMatrix100.size() > 0;
    }

    [[nodiscard]] bool integrationVecIsValid() const {
        return integrationVec.size() > 0;
    }
};

namespace {

    void filterXICPointsByAccuracyPPM(
        float mzVal,
        float ppmTol,
        XICPoints *xicPoints
        ) {

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTol);

        const auto terminatorLogic = [mzVal, massTol](const XICPoint &p) {
            return !(mzVal - massTol < p.mz && p.mz < mzVal + massTol);
        };

        const auto terminator = std::remove_if(
            xicPoints->begin(),
            xicPoints->end(),
            terminatorLogic
            );

        xicPoints->erase(terminator, xicPoints->end());
    }

    Err getXICs(
        const QVector<MS2Ion> &ms2Ions,
        float ppmTol,
        XICPeakManager *xicPeakManager,
        QVector<XICPoints> *xicPointsVec100,
        QVector<XICPoints> *xicPointsVec100Shadows,
        QVector<XICPoints> *xicPointsVec45,
        QVector<XICPoints> *xicPointsVec20
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;

        xicPointsVec100->resize(ms2Ions.size());

        for (int i = 0; i < ms2Ions.size(); i++) {

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

            if (xicPoints.empty()) {
                xicPointsVec100->push_back({});
                continue;
            }

            xicPointsVec100->push_back(xicPoints);

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
                &xicPoints
                );
            xicPointsVec45->push_back(xicPoints);

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
                &xicPoints
                );
            xicPointsVec20->push_back(xicPoints);

            XICPoints xicPointsShadows;
            const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
            e = xicPeakManager->getXIC(ms2Ion.mz - isotopeDistanceThomsons, &xicPointsShadows); ree;
            xicPointsVec100Shadows->push_back(xicPointsShadows);
        }

        ERR_RETURN
    }

    FrameIndex findFrameIndexMaxXICPointsVec (const QVector<XICPoints> &xicPointsVec100) {

        FrameIndex frameIndexMax = -1;

        for (const XICPoints &xicPoints : xicPointsVec100) {

            if (xicPoints.empty()) {
                continue;
            }

            const FrameIndex frameIndexMaxXICPoints = std::max_element(
                xicPoints.begin(),
                xicPoints.end(),
                [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
                )->scanNumber;
            frameIndexMax = std::max(frameIndexMax, frameIndexMaxXICPoints);
        }

        return frameIndexMax;
    }

    Err buildEigenMatrix(
        const QVector<XICPoints> &xicPointsVec,
        const Eigen::VectorX<float> &kernelMs2,
        FrameIndex frameIndexMax,
        Eigen::MatrixX<float> *matIntensity,
        Eigen::MatrixX<float> *matMz
        ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(xicPointsVec); ree;

        const FrameIndex frameIndexBuffer = 2;

        matIntensity->resize(frameIndexMax, xicPointsVec.size() + frameIndexBuffer);
        matIntensity->setZero();

        matMz->resize(frameIndexMax, xicPointsVec.size() + frameIndexBuffer);
        matMz->setZero();

        for (int col = 0; col < xicPointsVec.size(); col++) {

            const XICPoints &xicPointsCol = xicPointsVec.at(col);
            for (const XICPoint &p : xicPointsCol) {
                matIntensity->coeffRef(p.scanNumber, col) = p.intensity;
                matMz->coeffRef(p.scanNumber, col) = p.mz;
            }
        }

        *matIntensity = EigenKernelUtils::applyKernelToEachColumnInMatrix(*matIntensity, kernelMs2);

        ERR_RETURN
    }

    Err buildIntegrationVector(
        const MatriciesAndVecs &matriciesAndVecs,
        const Eigen::VectorX<float> &kernelIntegration,
        Eigen::VectorX<float> *integrationVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;
        constexpr float minPeakCount = 3.9;

        Eigen::MatrixX<float> matCount = matriciesAndVecs.intensityMatrix100;
        matCount = (matCount.array() > intensityThresholdVal).select(countValue, matCount);
        EigenUtils::thresholdMatrix(0.0f, &matCount);

        Eigen::VectorX<float> integrationVecLocal = matCount.rowwise().sum();
        EigenUtils::thresholdVector(minPeakCount, &integrationVecLocal);

        *integrationVec = EigenKernelUtils::convolveVectorWithKernel(
            integrationVecLocal,
            kernelIntegration
            );

        ERR_RETURN
    }

    Err initMatricesdAndVecs(
        const QVector<MS2Ion> &ms2Ions,
        const Eigen::VectorX<float> &kernelMs2,
        const Eigen::VectorX<float> &kernelIntegration,
        int topNMS2Ions,
        float ppmTol,
        XICPeakManager *xicPeakManager,
        MatriciesAndVecs *matriciesAndVecs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;

        QVector<MS2Ion> ms2IonsResized = ms2Ions;
        ms2IonsResized.resize(topNMS2Ions);

        QVector<XICPoints> xicPointsVec100;
        QVector<XICPoints> xicPointsVec100Shadow;
        QVector<XICPoints> xicPointsVec45;
        QVector<XICPoints> xicPointsVec20;
        e = getXICs(
            ms2Ions,
            ppmTol,
            xicPeakManager,
            &xicPointsVec100,
            &xicPointsVec100Shadow,
            &xicPointsVec45,
            &xicPointsVec20
            ); ree;

        const FrameIndex frameIndexMax = findFrameIndexMaxXICPointsVec(xicPointsVec100);

        e = buildEigenMatrix(
            xicPointsVec100,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix100,
            &matriciesAndVecs->mzMatrix100
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec45,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix45,
            &matriciesAndVecs->mzMatrix45
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec20,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix20,
            &matriciesAndVecs->mzMatrix20
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec100Shadow,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix100Shadow,
            &matriciesAndVecs->mzMatrix100Shadow
            ); ree;

        e = buildIntegrationVector(
            *matriciesAndVecs,
            kernelIntegration,
            &matriciesAndVecs->integrationVec
            ); ree;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::calculateScores(
    const QVector<MS2Ion> &ms2Ions,
    TargetDecoyCandidatePair* tdcp,
    CandidateScores *candidateScores
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2Ions); ree;

    MatriciesAndVecs matriciesAndVecs;
    e = initMatricesdAndVecs(
        ms2Ions,
        d_ptr->m_kernelMs2,
        d_ptr->m_kernelIntegration,
        m_topNMS2Ions,
        static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM),
        m_xicPeakManager,
        &matriciesAndVecs
        ); ree;

    if (m_msCalibratomatic.isInit()) {
        //TODO write code to truncate vector based on predicted scan time/frame index
    }

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        matriciesAndVecs.integrationVec,
        m_pythiaParameters.stopThresholdFraction,
        &peakIntegrationsVsIntensities
        ); ree;

    e = processIntegrationVectorPeakIntegrations(
        matriciesAndVecs,
        peakIntegrationsVsIntensities
        ); ree;

    ERR_RETURN
}

namespace {

    QPair<PeakIntegrationIndexes, Intensity> correctPeakIntegrationForSingleRow(
        const QPair<PeakIntegrationIndexes, Intensity> &pii,
        const Eigen::VectorX<float> &integrationVec
        ) {
        QPair<PeakIntegrationIndexes, Intensity> piiWorking = pii;
        if (piiWorking.first.first == piiWorking.first.second) {
            piiWorking.first.first = std::max(0, piiWorking.first.first - 1);
            piiWorking.first.second = std::min(
                static_cast<int>(integrationVec.size()) - 1,
                piiWorking.first.second + 1
                );
        }
        return piiWorking;
    }

}//namespace
Err CandidateScorertron::processIntegrationVectorPeakIntegrations(
    const MatriciesAndVecs &matriciesAndVecs,
    const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
    e = ErrorUtils::isTrue(matriciesAndVecs.integrationVecIsValid()); ree;

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensityResized = peakIntegrationsVsIntensity;

    constexpr int topNIntegrations = 10; //TODO make this settable
    peakIntegrationsVsIntensityResized.resize(topNIntegrations);

    for (const QPair<PeakIntegrationIndexes, Intensity> &pii : peakIntegrationsVsIntensity) {

        const QPair<PeakIntegrationIndexes, Intensity> piiWorking = correctPeakIntegrationForSingleRow(
            pii,
            matriciesAndVecs.integrationVec
            );

        Eigen::MatrixX<float> mat = m_pythiaParameters.subtractShadows
                                  ? matriciesAndVecs.intensityMatrix100 - matriciesAndVecs.intensityMatrix100Shadow
                                  : matriciesAndVecs.intensityMatrix100;
        EigenUtils::thresholdMatrix(0.0f, &mat);

        Eigen::MatrixX<float> matBlock = mat.block(
              piiWorking.first.first,
              0,
              piiWorking.first.second - piiWorking.first.first + 1,
              mat.cols()
              ).eval();


    }





    ERR_RETURN
}



