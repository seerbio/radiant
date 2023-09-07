//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "IRTPredictron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "ParallelUtils.h"
#include "ScanTimeFromIRtMapper.h"
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

    m_params = params;
    m_fragPredsTopN = peptideStringWithModsVsCandidatePeptide;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = FragLibReader::generateFragmentFrequencies(
            m_fragPredsTopN,
            m_params.ms2ExtractionWidthPPM,
            &m_fragmentFrequencies
    ); ree

    e = m_msFrame.init(
            scanNumberVsScanPoints,
            scanNumberVsScanTime
            ); ree;


//NOTE: Turn off deisotoping in PythiaDIAWorkflow.cpp if using here.
//#define DEISOTOPE
#ifdef DEISOTOPE
    e = m_msFrame.deisotopeMsFrame(m_params.ms2ExtractionWidthPPM); ree;
#endif

    e = m_msFrameMS1.init(
            scanNumberVsScanPointsMS1,
            scanNumberVsScanTime
            ); ree;

    ERR_RETURN
}

Err MsFrameScoretron::init(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const PythiaParameters &params,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QString &iRTRecalibrationFilePath
        ) {

    ERR_INIT

    e = init(
            uniqueMsInfoScanKey,
            params,
            scanNumberVsScanPoints,
            scanNumberVsScanPointsMS1,
            peptideStringWithModsVsCandidatePeptide,
            scanNumberVsScanTime
            ); ree;

    qDebug() << "updating rt vals from iRT";

    ScanTimeFromIRtMapper scanTimeFromIRtMapper;
    e = scanTimeFromIRtMapper.init(iRTRecalibrationFilePath); ree;

    for (auto it = peptideStringWithModsVsCandidatePeptide.begin(); it != peptideStringWithModsVsCandidatePeptide.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const double iRT = it.value().iRt;

        double predictedScanTime;
        e = scanTimeFromIRtMapper.predictScanTime(iRT, &predictedScanTime); ree;

        m_fragPredsPredictedScanTime.insert(peptideStringWithMods, predictedScanTime);
    }

    ERR_RETURN
}

Err MsFrameScoretron::scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates) {

    ERR_INIT

    QMap<PeptideStringWithMods, CandidatePeptide> peptideStringWithModsVsCandidatePeptideDecoys;
    e = buildPeptideStringWithModsVsCandidatePeptideDecoys(
            &peptideStringWithModsVsCandidatePeptideDecoys
            ); ree;

    QMap<MzHashed, XICPoints> mzHashedVsXICPoints;
    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
    e = buildMS2Peaks(
            m_fragPredsTopN,
            &mzHashedVsXICPoints,
            &mzHashedVsIonPresence
            ); ree;

    e = buildMS2Peaks(
            peptideStringWithModsVsCandidatePeptideDecoys,
            &mzHashedVsXICPoints,
            &mzHashedVsIonPresence
    ); ree;

    if (m_fragPredsPredictedScanTime.isEmpty()) {
        e = m_candidateProcessertron.init(
                m_params,
                mzHashedVsXICPoints,
                mzHashedVsIonPresence,
                m_msFrame,
                m_msFrameMS1,
                m_scanNumberVsScanTime,
                m_fragmentFrequencies,
                m_uniqueMsInfoScanKey
        ); ree;
    }
    else {
        e = m_candidateProcessertron.init(
                m_params,
                mzHashedVsXICPoints,
                mzHashedVsIonPresence,
                m_msFrame,
                m_msFrameMS1,
                m_scanNumberVsScanTime,
                m_fragmentFrequencies,
                m_uniqueMsInfoScanKey,
                m_fragPredsPredictedScanTime
        ); ree;
    }


    for (const CandidatePeptide &candidatePeptide : m_fragPredsTopN) {

        ScoredCandidate scoredCandidate;
        e = m_candidateProcessertron.processCandidateTarget(
                candidatePeptide,
                &scoredCandidate
                ); ree;

        if (scoredCandidate.cosineSimSum < 0.0) {
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

namespace {

    Err mutateCandidatePeptideTarget(
            const CandidatePeptide &candidatePeptideTarget,
            CandidatePeptide *candidatePeptideDecoy
    ) {

        ERR_INIT

        const QMap<QChar, double> diannMutateAminoAcidTo = AminoAcids::diannMutateAminoAcidTo();

        const QString &seq = candidatePeptideTarget.peptideStringWithMods;

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = seq.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        QVector<MS2Ion> ms2IonDecoys;
        for (const MS2Ion &ms2Ion : candidatePeptideTarget.ms2Ions) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }
                }
            }

            else if (ionLableInfo.second.contains('y')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }
                }
            }

            else {
                qDebug() << "Non b/y ion" << ionLableInfo;
            }

            ms2IonDecoys.push_back(ms2IonDecoy);
        }

        *candidatePeptideDecoy = candidatePeptideTarget;
        candidatePeptideDecoy->isDecoy = true;
        candidatePeptideDecoy->mass += nTermDeltaMass + cTermDeltaMass;
        candidatePeptideDecoy->ms2Ions = ms2IonDecoys;

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::buildPeptideStringWithModsVsCandidatePeptideDecoys(
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptideDecoys
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPredsTopN); ree;

    for (const CandidatePeptide &candidatePeptideTarget : m_fragPredsTopN) {

        CandidatePeptide candidatePeptideDecoy;
        e = mutateCandidatePeptideTarget(
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

    Err buildMzHashedVsMzIon(
            const QMap<PeptideStringWithMods, CandidatePeptide> &fragPredsTopN,
            QMap<MzHashed, MZION> *mzHashedVsMzIon
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(fragPredsTopN); ree;
        mzHashedVsMzIon->clear();

        for (auto it = fragPredsTopN.begin(); it != fragPredsTopN.end(); it++) {

            const QVector<MS2Ion> &ms2IonsTopN = it.value().ms2Ions;

            for (const MS2Ion &ms2Ion: ms2IonsTopN) {
                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                mzHashedVsMzIon->insert(mzHashed, ms2Ion.mz);
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
            vecPresence = EigenUtils::setNANToZero(vecPresence);
            mzHashedVsIonPresence->insert(mzHashed, EigenUtils::convertEigenVectorToQVector(vecPresence));
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::buildMS2Peaks(
        const QMap<PeptideStringWithMods, CandidatePeptide> &candidatePeptides,
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints,
        QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidatePeptides); ree;

    QMap<MzHashed, MZION> mzHashedVsMzIon;
    e = buildMzHashedVsMzIon(candidatePeptides, &mzHashedVsMzIon); ree;

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM,
            mzHashedVsXICPoints
            ); ree;

    QMap<MzHashed, Eigen::VectorX<double>> mzHashedVsXICPointsNormalized;
    e = buildMzHashedVsXICPointsNormalized(
            *mzHashedVsXICPoints,
            m_msFrame.scanCount(),
            &mzHashedVsXICPointsNormalized
            ); ree;

    e = buildMzHashedVsIonPresence(
            mzHashedVsXICPointsNormalized,
            mzHashedVsIonPresence
            );ree

    ERR_RETURN
}


