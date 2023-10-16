//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "IRTPredictron.h"
#include "MsCalibratomatic.h"
#include "ParallelUtils.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"


namespace{

    PeakIntegratomaticParameters buildPeakIntegratomaticParams(const PythiaParameters  &pythiaParameters) {

        PeakIntegratomaticParameters params;
        params.smoothCount = pythiaParameters.smoothCount;
        params.filterLength = pythiaParameters.filterLength;
        params.sigma = pythiaParameters.sigma;
        params.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        return params;
    }

}//namespace
Err MsFrameScoretron::init(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const PythiaParameters &params,
        int topNMS2Ions,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanPointsMS1); ree;
    e = ErrorUtils::isNotEmpty(peptideStringWithModsVsCandidatePeptide); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isTrue(topNMS2Ions >= 6); ree;

    m_params = params;
    m_fragPredsTopN = peptideStringWithModsVsCandidatePeptide;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_scanNumberVsScanTime = scanNumberVsScanTime;
    m_topNMS2Ions = topNMS2Ions;

    e = m_msFrame.init(
            scanNumberVsScanPoints,
            scanNumberVsScanTime
            ); ree;

    e = m_msFrameMS1.init(
            scanNumberVsScanPointsMS1,
            scanNumberVsScanTime
            ); ree;

    ERR_RETURN
}

namespace {

    Err buildFragPredsPredictedScanTime(
            const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
            MsCalibratomatic *msCalibratomatic,
            QMap<PeptideStringWithMods, ScanTime> *fragPredsPredictedScanTime
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsCandidatePeptide); ree;
        e = ErrorUtils::isTrue(msCalibratomatic->isInit()); ree;

        for (auto it = peptideStringWithModsVsCandidatePeptide.begin(); it != peptideStringWithModsVsCandidatePeptide.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const double iRT = it.value().iRt;

            double predictedScanTime;
            e = msCalibratomatic->predictScanTime(iRT, &predictedScanTime); ree;

            fragPredsPredictedScanTime->insert(peptideStringWithMods, predictedScanTime);
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const PythiaParameters &params,
        int topNMS2Ions,
        const QMap<ScanNumber, ScanPoints> &_scanNumberVsScanPoints,
        const QMap<ScanNumber, ScanPoints> &_scanNumberVsScanPointsMS1,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const MsCalibratomatic &msCalibratomatic
        ) {

    ERR_INIT

    m_msCalibratomatic = msCalibratomatic;

    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints = _scanNumberVsScanPoints;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1 = _scanNumberVsScanPointsMS1;
    if (m_msCalibratomatic.isInit()) {

        e = m_msCalibratomatic.recalibrateScanPoints(
                _scanNumberVsScanPoints,
                &scanNumberVsScanPoints
                ); ree;

        e = m_msCalibratomatic.recalibrateScanPoints(
                _scanNumberVsScanPointsMS1,
                &scanNumberVsScanPointsMS1
        ); ree;

        e = buildFragPredsPredictedScanTime(
                peptideStringWithModsVsCandidatePeptide,
                &m_msCalibratomatic,
                &m_fragPredsPredictedScanTime
        ); ree;
    }

    e = init(
            uniqueMsInfoScanKey,
            params,
            topNMS2Ions,
            scanNumberVsScanPoints,
            scanNumberVsScanPointsMS1,
            peptideStringWithModsVsCandidatePeptide,
            scanNumberVsScanTime
    ); ree;


    ERR_RETURN
}

Err MsFrameScoretron::scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates) {

    ERR_INIT

    QMap<PeptideStringWithMods, CandidatePeptide> peptideStringWithModsVsCandidatePeptideDecoys;
    e = buildPeptideStringWithModsVsCandidatePeptideDecoys(
            &peptideStringWithModsVsCandidatePeptideDecoys
            ); ree;

    QMap<MzHashed, XICPoints> mzHashedVsXICPoints100;
    QMap<MzHashed, XICPoints> mzHashedVsXICPoints100Shadows;
    QMap<MzHashed, XICPoints> mzHashedVsXICPoints45;
    QMap<MzHashed, XICPoints> mzHashedVsXICPoints20;
    QMap<MzHashed, XICPoints> mzHashedVsXICPointsB2B3;
    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
    e = buildMS2Peaks(
            m_fragPredsTopN,
            &mzHashedVsXICPoints100,
            &mzHashedVsXICPoints100Shadows,
            &mzHashedVsXICPoints45,
            &mzHashedVsXICPoints20,
            &mzHashedVsXICPointsB2B3,
            &mzHashedVsIonPresence
            ); ree;

    e = buildMS2Peaks(
            peptideStringWithModsVsCandidatePeptideDecoys,
            &mzHashedVsXICPoints100,
            &mzHashedVsXICPoints100Shadows,
            &mzHashedVsXICPoints45,
            &mzHashedVsXICPoints20,
            &mzHashedVsXICPointsB2B3,
            &mzHashedVsIonPresence
    ); ree;

    if (m_fragPredsPredictedScanTime.isEmpty()) {
        e = m_candidateProcessertron.init(
                m_params,
                m_topNMS2Ions,
                mzHashedVsXICPoints100,
                mzHashedVsXICPoints100Shadows,
                mzHashedVsXICPoints45,
                mzHashedVsXICPoints20,
                mzHashedVsXICPointsB2B3,
                mzHashedVsIonPresence,
                m_msFrame,
                m_msFrameMS1,
                m_scanNumberVsScanTime,
                m_uniqueMsInfoScanKey
        ); ree;
    }
    else {
        e = m_candidateProcessertron.init(
                m_params,
                m_topNMS2Ions,
                mzHashedVsXICPoints100,
                mzHashedVsXICPoints100Shadows,
                mzHashedVsXICPoints45,
                mzHashedVsXICPoints20,
                mzHashedVsXICPointsB2B3,
                mzHashedVsIonPresence,
                m_msFrame,
                m_msFrameMS1,
                m_scanNumberVsScanTime,
                m_uniqueMsInfoScanKey,
                m_fragPredsPredictedScanTime
        ); ree;
    }

    const double cosineSimSumMin = 0.0;
    for (const CandidatePeptide &candidatePeptide : m_fragPredsTopN) {

        ScoredCandidate scoredCandidate;
        e = m_candidateProcessertron.processCandidateTarget(
                candidatePeptide,
                &scoredCandidate
                ); ree;

        if (scoredCandidate.cosineSimSum100 < cosineSimSumMin) {
            continue;
        }

        scoredCandidates->push_back(scoredCandidate);
    }

    QMap<PeptideStringWithMods, ScoredCandidate> scoredCandidateDecoys;
    for (const ScoredCandidate &scoredCandidateTarget : *scoredCandidates) {

        e = ErrorUtils::isTrue(
                peptideStringWithModsVsCandidatePeptideDecoys.contains(scoredCandidateTarget.peptideStringWithMods)
                );ree

        ScoredCandidate scoredCandidateDecoy;
        e = m_candidateProcessertron.processCandidateDecoy(
                peptideStringWithModsVsCandidatePeptideDecoys.value(scoredCandidateTarget.peptideStringWithMods),
                scoredCandidateTarget.scanTime,
                &scoredCandidateDecoy
        ); ree;

        scoredCandidateDecoys.insert(scoredCandidateDecoy.peptideStringWithMods, scoredCandidateDecoy);
    }

    scoredCandidates->append(scoredCandidateDecoys.values().toVector());

    ERR_RETURN
}

Err MsFrameScoretron::buildPeptideStringWithModsVsCandidatePeptideDecoys(
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptideDecoys
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPredsTopN); ree;

    for (const CandidatePeptide &candidatePeptideTarget : m_fragPredsTopN) {

        CandidatePeptide candidatePeptideDecoy;
        e = FragLibReader::mutateCandidatePeptideTarget(
                candidatePeptideTarget,
                &candidatePeptideDecoy
        ); ree;

        peptideStringWithModsVsCandidatePeptideDecoys->insert(
                candidatePeptideTarget.peptideStringWithMods,
                candidatePeptideDecoy
                );
    }

    ERR_RETURN
}

//TODO make buildMS2Peaks() into its own class
namespace {

    enum class IonSelector {
        MS2Ions,
        MS2ShadowIons,
        B2B3Ions,
        Y2Y3Ions
    };

    Err buildMzHashedVsMzIon(
            const QMap<PeptideStringWithMods, CandidatePeptide> &fragPredsTopN,
            IonSelector ionSelector,
            QMap<MzHashed, MZION> *mzHashedVsMzIon
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(fragPredsTopN); ree;
        mzHashedVsMzIon->clear();

        for (auto it = fragPredsTopN.begin(); it != fragPredsTopN.end(); it++) {

            if (ionSelector == IonSelector::Y2Y3Ions) {
                //removed code that was. Not collecting data for Y2Y3 anymore
            }
            else if (ionSelector == IonSelector::B2B3Ions) {
                const QVector<MZION> &mzIons = it.value().ms2IonMzB2B3;
                for (const MZION &mz: mzIons) {
                    const MzHashed mzHashed = MathUtils::hashDecimal(mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                    mzHashedVsMzIon->insert(mzHashed, mz);
                }
            }
            else if (ionSelector == IonSelector::MS2ShadowIons) {

                const QVector<MS2Ion> &ms2IonsTopN = it.value().ms2Ions;
                const int maxShadowCount = 6;

                int counter = 0;
                for (const MS2Ion &ms2Ion: ms2IonsTopN) {

                    if (counter++ >= maxShadowCount) {
                        break;
                    }

                    e = ErrorUtils::isTrue(ms2Ion.charge > 0); ree;
                    const double isotopeChargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
                    const double mzIsotopeShadow = ms2Ion.mz - isotopeChargeDistance;

                    const MzHashed mzHashed = MathUtils::hashDecimal(mzIsotopeShadow, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                    mzHashedVsMzIon->insert(mzHashed, mzIsotopeShadow);
                }
            }
            else {
                const QVector<MS2Ion> &ms2IonsTopN = it.value().ms2Ions;
                for (const MS2Ion &ms2Ion: ms2IonsTopN) {
                    const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                    mzHashedVsMzIon->insert(mzHashed, ms2Ion.mz);
                }
            }
        }

        ERR_RETURN
    }

    Err buildMzHashedVsXICPoints(
            const QMap<MzHashed, MZION> &mzHashedVsMzIon,
            const MsFrame &msFrame,
            double ppmTol,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints
            ) {

        ERR_INIT

        TurboXIC turboXic;
        e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

        double scanNumberMin;
        double scanNumberMax;
        double mzMinRTree;
        double mzMaxRTree;
        e = turboXic.getRTreeLimits(
                &scanNumberMin,
                &scanNumberMax,
                &mzMinRTree,
                &mzMaxRTree
        ); ree;

        for (auto it = mzHashedVsMzIon.begin(); it != mzHashedVsMzIon.end(); it++) {

            const MzHashed mzHashed = it.key();
            if (mzHashedVsXICPoints->contains(mzHashed)) {
                continue;
            }

            const MZION mz = it.value();
            const double mzTol = MathUtils::calculatePPM(mz, ppmTol);

            const double mzMin = mz - mzTol;
            const double mzMax = mz + mzTol;

            const XICPoints xicPoints = turboXic.extractPointsXIC(
                    mzMin,
                    mzMax,
                    static_cast<int>(scanNumberMin),
                    static_cast<int>(scanNumberMax)
            );

            mzHashedVsXICPoints->insert(mzHashed, xicPoints);

        }

        ERR_RETURN
    }

    Err buildMzHashedVsXICPointsNormalized(
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            int scanCount,
            QMap<MzHashed, Eigen::VectorX<double>> *mzHashedVsXICPointsNormalized
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;

        for (auto it = mzHashedVsXICPoints.begin(); it != mzHashedVsXICPoints.end(); it++) {

            const MzHashed mzHashed = it.key();
            const XICPoints &xicPoints = it.value();

            if (xicPoints.scanNumbersVsIntensityVals.isEmpty()) {
                mzHashedVsXICPointsNormalized->insert(mzHashed, {});
                continue;
            }

            Eigen::VectorX<double> vecNormalized = EigenUtils::convertQMapToEigenVector(
                    xicPoints.scanNumbersVsIntensityVals,
                    scanCount
            );
            const double denom = vecNormalized.maxCoeff();
            vecNormalized = vecNormalized.array() / denom;

            mzHashedVsXICPointsNormalized->insert(mzHashed, vecNormalized);
        }

        ERR_RETURN
    }

    Err buildMzHashedVsIonPresence(
            const QMap<MzHashed, Eigen::VectorX<double>> &mzHashedVsXICPointsNormalized,
            QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPointsNormalized); ree;

        for (auto it = mzHashedVsXICPointsNormalized.begin(); it != mzHashedVsXICPointsNormalized.end(); it++) {

            const MzHashed mzHashed = it.key();
            const Eigen::VectorX<double> &normVec = it.value();

            Eigen::VectorX<double> vecPresence = normVec.array() / normVec.array();
            EigenUtils::replaceNaN(0.0, &vecPresence);
            mzHashedVsIonPresence->insert(mzHashed, EigenUtils::convertEigenVectorToQVector(vecPresence));
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::buildMS2Peaks(
        const QMap<PeptideStringWithMods, CandidatePeptide> &candidatePeptides,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints100,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints100Shadows,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints45,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints20,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPointsB2B3,
        QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidatePeptides); ree;

    QMap<MzHashed, MZION> mzHashedVsMzIon;
    e = buildMzHashedVsMzIon(
            candidatePeptides,
            IonSelector::MS2Ions,
            &mzHashedVsMzIon
            ); ree;

    QMap<MzHashed, MZION> mzHashedVsMzIonB2B3;
    e = buildMzHashedVsMzIon(
            candidatePeptides,
            IonSelector::B2B3Ions,
            &mzHashedVsMzIonB2B3
    ); ree;

    QMap<MzHashed, MZION> mzShadowsHashedVsMzIon;
    e = buildMzHashedVsMzIon(
            candidatePeptides,
            IonSelector::MS2ShadowIons,
            &mzShadowsHashedVsMzIon
    ); ree;

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM,
            mzHashedVsXICPoints100
            ); ree;

    e = buildMzHashedVsXICPoints(
            mzShadowsHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM,
            mzHashedVsXICPoints100Shadows
    ); ree;

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIonB2B3,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM,
            mzHashedVsXICPointsB2B3
    ); ree;

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
            mzHashedVsXICPoints45
    ); ree;

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
            mzHashedVsXICPoints20
    ); ree;

    QMap<MzHashed, Eigen::VectorX<double>> mzHashedVsXICPoints100Normalized;
    e = buildMzHashedVsXICPointsNormalized(
            *mzHashedVsXICPoints100,
            m_msFrame.scanCount(),
            &mzHashedVsXICPoints100Normalized
            ); ree;

    e = buildMzHashedVsIonPresence(
            mzHashedVsXICPoints100Normalized,
            mzHashedVsIonPresence
            );ree

    ERR_RETURN
}
