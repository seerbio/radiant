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
#include "MsCalibratomatic.h"
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
            int topNMS2Ions,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
            const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const MsCalibratomatic &msCalibratomatic
    );

    Err init(
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const PythiaParameters &params,
            int topNMS2Ions,
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
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints100,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints100Shadows,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints45,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints20,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPointsB2B3,
            QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
            );


private:

    QMap<PeptideStringWithMods, CandidatePeptide> m_fragPredsTopN;
    QMap<PeptideStringWithMods, ScanTime> m_fragPredsPredictedScanTime;

    PythiaParameters m_params;
    int m_topNMS2Ions;
    QString m_msDataFilePath;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;
    MsCalibratomatic m_msCalibratomatic;
    MsFrame m_msFrame;
    MsFrame m_msFrameMS1;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    CandidateProcessertron m_candidateProcessertron;

    MS2ChargeDeconvolvotron m_ms2ChargeDeconvolvotron;
};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
