//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "CalibrationReader.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FrameExtractsReader.h"
#include "MsFrameScoretronProcessormatic.h"
#include "NearestNeighborsSearch.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>


Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
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
            &m_fragPredsIsDecoy,
            &m_fragPredsMass,
            &m_fragPredsIRT
    ); ree

    e = loadFragPredsTopN(m_params.topNMs2Ions); ree

    e = initMsFrame(
        msDataFilePath,
        uniqueMsInfoScanKey,
        mzTargetStartStop
        ); ree

    m_mzStartStopMean = (mzTargetStartStop.second + mzTargetStartStop.first) / 2.0;

//    const QString &chargeModelFilePath
//            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");
//
//    const QString &monoModelFilePath
//            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");
//
//    e = m_ms2ChargeDeconvolvotron.init(
//            chargeModelFilePath,
//            monoModelFilePath,
//            m_params.ms2ExtractionWidthPPM
//            ); ree

    qDebug() << "TargetKey" << uniqueMsInfoScanKey;
    qDebug() << "Scan Count" << m_msFrame.scanCount();
    qDebug() << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}

namespace {

    Err buildNearestNeighborsIRTData(
            const QString &iRTRecalibrationFilePath,
            QVector<QPair<double, Coors>> *nnInputData
    ) {

        ERR_INIT

        e = ErrorUtils::fileExists(iRTRecalibrationFilePath); ree

        QVector<IRTReCalibrationRow> iRTReCalibrationReaderRows;
        e  = CSVReader::read(
                iRTRecalibrationFilePath,
                &iRTReCalibrationReaderRows
        ); ree;

        for (const IRTReCalibrationRow &row : iRTReCalibrationReaderRows) {
            nnInputData->push_back({row.scanTime, {row.iRT, 0.0}});
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const QString &iRTRecalibrationFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
        ) {
    ERR_INIT

    e = init(
            params,
            msDataFilePath,
            fragLibFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop
            ); ree;

    qDebug() << "updating rt vals from iRT";

    QVector<QPair<double, Coors>> nnInputData;
    e = buildNearestNeighborsIRTData(
            iRTRecalibrationFilePath,
            &nnInputData
    );

    NearestNeighborsSearch nearestNeighborsSearch;
    e = nearestNeighborsSearch.init(nnInputData); ree

    const int kNearestPoints = 10;

    for (auto it = m_fragPredsIRT.begin(); it != m_fragPredsIRT.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const double iRT = it.value();

        QVector<NNSearchResult> nnSearchResult;
        nearestNeighborsSearch.kNearestNeighborsSearch(
                {{iRT, 0.0}},
                kNearestPoints,
                &nnSearchResult
        );

        m_fragPredsPredictedScanTime.insert(peptideStringWithMods, nnSearchResult.front().meanValues);
    }

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

    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
    );ree;

    ERR_RETURN
}

Err MsFrameScoretron::loadFragPredsTopN(int n) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPreds); ree;

    for (auto it = m_fragPreds.begin(); it != m_fragPreds.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const MS2IonsSeparated &ms2IonsSeparated = it.value();

        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.y2Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.b2Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.aIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yNH3Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yH2OIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bNH3Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bH2OIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.precursorIons.values().toVector());

        QVector<MS2Ion> *ms2Ions = &m_fragPredsTopN[peptideStringWithMods];
        FragLibReader::getTopNMostIntenseMs2Ions(n, ms2Ions);
    }

    ERR_RETURN
}

Err MsFrameScoretron::scoreFrameCandidates(QMap<FrameIndex , QVector<ScoredCandidate>> *frameIndexVsScoredCandidates) {

    ERR_INIT


//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    for(auto it = apexScanPoints.begin(); it != apexScanPoints.end(); it++) {
        const ScanNumber sn = m_msFrame.scanNumberFromFrameIndex(it.key());
        scanNumberVsScanPoints.insert(sn, it.value());
    }

    e = MsFrame::writeFrameScans(scanNumberVsScanPoints, "frame.prq"); ree
#endif

    ERR_RETURN
}

//namespace{
//
//    Err buildScanFragPreds(
//            const QVector<ScoredCandidate> &scoredCandidatesTargets,
//            QMap<PeptideStringWithMods, QVector<MS2Ion>> *fragPredsFlattened,
//            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanFragPredsFlattened
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(*fragPredsFlattened); ree;
//
//        for (const ScoredCandidate &sc : scoredCandidatesTargets) {
//
//            e = ErrorUtils::isTrue(fragPredsFlattened->contains(sc.peptideStringWithMods)); ree;
//
//            const QVector<MS2Ion> &pred = fragPredsFlattened->value(sc.peptideStringWithMods);
//
//            scanFragPredsFlattened->insert(sc.peptideStringWithMods, pred);
//        }
//
//        ERR_RETURN
//    }
//
//    void removeNegativeDiscScores(QMap<PeptideStringWithMods, TandemDeconvolverResult> *pepSeqVsWeight) {
//
//        QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeightFiltered;
//        for (auto it = pepSeqVsWeight->begin(); it != pepSeqVsWeight->end(); it++) {
//
//            const PeptideStringWithMods &peptideStringWithMods = it.key();
//            const TandemDeconvolverResult &res = it.value();
//
//            if (res.discScore < 0 || res.pVal > 0.05) {
//                continue;
//            }
//
//            pepSeqVsWeightFiltered.insert(peptideStringWithMods, res);
//        }
//
//        *pepSeqVsWeight = pepSeqVsWeightFiltered;
//    }
//
//}//namespace
//Err MsFrameScoretron::deconvolveScan(
//        const QVector<ScoredCandidate> &scoredCandidatesTargets,
//        const ScanPoints &scanPointsDeisotoped,
//        QVector<ScoredCandidate> *scoredCandidatesTargetsFiltered
//        ) {
//
//    ERR_INIT
//
//    QMap<PeptideStringWithMods, QVector<MS2Ion>> scanFragPredsTopN;
//    e = buildScanFragPreds(
//            scoredCandidatesTargets,
//            &m_fragPredsTopN,
//            &scanFragPredsTopN
//            ); ree
//
//    const int precision = 2;
//    const double mzMax = 1400.0;
//    const int iters = 4;
//    const double stopTol = 1e-8;
//    const double pVal = 0.05;
//
//    if (scanFragPredsTopN.isEmpty()) {
//        ERR_RETURN
//    }
//
//    TandemSpectraDeconvolvotron decon;
//    e = decon.init(
//            precision,
//            mzMax,
//            iters,
//            stopTol,
//            pVal
//            ); ree
//
//    QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeight;
//    e = decon.deconvolveTandemSpectra(scanPointsDeisotoped, scanFragPredsTopN, &pepSeqVsWeight); ree
//
//    removeNegativeDiscScores(&pepSeqVsWeight);
//
//    for (const ScoredCandidate &sc : scoredCandidatesTargets) {
//
//        if (!pepSeqVsWeight.contains(sc.peptideStringWithMods)) {
//            continue;
//        }
//
//        const TandemDeconvolverResult &tdr = pepSeqVsWeight.value(sc.peptideStringWithMods);
//
//        ScoredCandidate scNew = sc;
//
//        scNew.discScore = tdr.discScore;
//        scNew.pVal = tdr.pVal;
//        scNew.frameError = tdr.frameError;
//        scNew.tTestVal = tdr.tTestVal;
//
//        scoredCandidatesTargetsFiltered->push_back(scNew);
//    }
//
//    ERR_RETURN
//}
