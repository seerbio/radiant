//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "CandidateProcessertron.h"
#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsFrame.h"
#include "ParquetReader.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"


using namespace Error;

class ScoredCandidate;
class TandemDeconvolverResult;
class XICPoints;

class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    friend class MissingPeptideManualTroubleshooter;
    friend class MsFrameScoretronProcessormaticTests;
    friend class MsFrameScoretronTests;

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    Err init(
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const PythiaParameters &params,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
            const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const QString &iRTRecalibrationFilePath
    );

    Err init(
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const PythiaParameters &params,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
            const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsMS2Ions,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
            );

    Err scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates);


private:

    Err buildPeptideStringWithModsVsCandidatePeptideDecoys(
            QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptideDecoys
            );

    Err buildMS2Peaks(
            const QMap<PeptideStringWithMods, CandidatePeptide> &candidatePeptides,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints,
            QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
            );


private:

    QMap<PeptideStringWithMods, CandidatePeptide> m_fragPredsTopN;
    QMap<MzHashed, FrequencyPercent> m_fragmentFrequencies;
    QMap<PeptideStringWithMods, ScanTime> m_fragPredsPredictedScanTime;

    PythiaParameters m_params;
    QString m_msDataFilePath;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;
    MsFrame m_msFrame;
    MsFrame m_msFrameMS1;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    CandidateProcessertron m_candidateProcessertron;

    MS2ChargeDeconvolvotron m_ms2ChargeDeconvolvotron;
};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
