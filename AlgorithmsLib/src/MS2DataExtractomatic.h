//
// Created by anichols on 9/30/23.
//

#ifndef PYTHIADIACPP_MS2DATAEXTRACTOMATIC_H
#define PYTHIADIACPP_MS2DATAEXTRACTOMATIC_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsReaderPointerAcc.h"
#include "ParquetReader.h"
#include "PythiaParameterReader.h"

using namespace Error;

class CandidatePeptide;
class ScoredCandidate;

class ALGORITHMSLIB_EXPORTS MS2DataExtractomatic {

public:

    MS2DataExtractomatic();
    ~MS2DataExtractomatic() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMs2Ions,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            MsReaderPointerAcc *msReaderPointerAcc
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMs2Ions,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            MsReaderPointerAcc *msReaderPointerAcc,
            const MsCalibratomatic &msCalibratomatic
    );

    Err extractMS2ForCandidates(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide,
            QVector<ScoredCandidate> *scoredCandidatesAll
            );

    static Err filterScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            bool filterDecoys,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
    );

    static Err outputFDRResults(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QMap<QString, int> *fdrVsCount
            );

private:

    Err extractTargetDecoyData(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptideCalibration,
            QVector<ScoredCandidate> *combinedResults
    );

    Err buildScoredCandidatesFDR(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            const QPair<double, double> &scanTimeMinMax,
            QVector<ScoredCandidate> *scoredCandidatesAll
    );

    Err fitWeightsLogic(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            const QPair<double, double> &scanTimeMinMax,
            const QVector<QVector<double>> &A,
            const QVector<double> &b,
            QVector<double> *weights,
            QVector<ScoredCandidate> *scoredCandidatesAll
    );

    static Err buildUniqueMsInfoScanKeyVsScanPoints(
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
            QMap<ScanNumber, ScanPoints > *scanNumberVsScanTimeMS1
    );

private:

    PythiaParameters m_pythiaParameters;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    int m_topNMs2Ions;

    MsCalibratomatic m_msCalibratomatic;

    bool m_useExtendedScores;
    bool m_useNeuralNetworkScores;


};


#endif //PYTHIADIACPP_MS2DATAEXTRACTOMATIC_H
