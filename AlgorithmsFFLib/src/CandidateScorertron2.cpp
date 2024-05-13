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

    e = d_ptr->init(); ree;

    ERR_RETURN
}

namespace {

    struct MatriciesAndVecs {

        Eigen::MatrixX<float> intensityMatrix100;
        Eigen::MatrixX<float> intensityMatrix100Shadow;
        Eigen::MatrixX<float> intensityMatrix45;
        Eigen::MatrixX<float> intensityMatrix20;


        Eigen::VectorX<float> integrationVec;
    };

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
        QVector<XICPoints> *xicPointsVec,
        QVector<XICPoints> *xicPointsVec45,
        QVector<XICPoints> *xicPointsVec20
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;

        xicPointsVec->resize(ms2Ions.size());

        for (int i = 0; i < ms2Ions.size(); i++) {

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

            if (xicPoints.empty()) {
                xicPointsVec->push_back({});
                continue;
            }

            xicPointsVec->push_back(xicPoints);

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
        }

        ERR_RETURN
    }

    Err initMatricesdAndVecs(
        const QVector<MS2Ion> &ms2Ions,
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

        QVector<XICPoints> xicPointsVec;
        QVector<XICPoints> xicPointsVec45;
        QVector<XICPoints> xicPointsVec20;
        e = getXICs(
            ms2Ions,
            ppmTol,
            xicPeakManager,
            &xicPointsVec,
            &xicPointsVec45,
            &xicPointsVec20
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
        m_topNMS2Ions,
        m_pythiaParameters.ms2ExtractionWidthPPM,
        m_xicPeakManager,
        &matriciesAndVecs
        ); ree;




    ERR_RETURN
}



