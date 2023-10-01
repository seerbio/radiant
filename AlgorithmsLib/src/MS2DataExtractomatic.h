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
#include "MsReaderParquet.h"
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
            MsReaderParquet *msReaderParquet
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMs2Ions,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            MsReaderParquet *msReaderParquet,
            const MsCalibratomatic &msCalibratomatic
    );

    Err extractMS2ForCandidates(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide,
            double fdrThreshold,
            QVector<ScoredCandidate> *scoredCandidatesAll,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
            );

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
            MsReaderParquet *msReaderParquet,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
            QMap<ScanNumber, ScanPoints > *scanNumberVsScanTimeMS1
    );

private:

    PythiaParameters m_pythiaParameters;
    MsReaderParquet *m_msReaderParquet;
    int m_topNMs2Ions;

    MsCalibratomatic m_msCalibratomatic;

    bool m_useExtendedScores;
    bool m_useNeuralNetworkScores;


};


#endif //PYTHIADIACPP_MS2DATAEXTRACTOMATIC_H
