//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;


class MsCalibratomatic;
class TargetDecoyPairParallelInput;
class TurboXIC;

class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairScoretron {

public:

    TargetDecoyCandidatePairScoretron();
    ~TargetDecoyCandidatePairScoretron() = default;

    /**
    * @brief Initializes the TargetDecoyCandidatePairScoretron.
    *
    * This method initializes the TargetDecoyCandidatePairScoretron with the provided parameters,
    * including Pythia parameters, MS1 scan information, MS reader accessor, DIA target frames,
    * and TurboXIC for MS1.
    *
    * @param pythiaParameters The Pythia parameters used for scoring.
    * @param scanNumberVsScanTimeMS1 Map of scan numbers to MS1 scan points for the scoring reference.
    * @param msReaderPointerAcc Pointer to the MS reader accessor for accessing MS data.
    * @param diaTargetFrames Map of DIA target frames for additional scoring context.
    * @param turboXICMS1 Pointer to TurboXIC for MS1, required for scoring calculations.
    * @return Error code indicating success or failure during initialization.
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanTimeMS1,
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
            TurboXIC *turboXICMS1
            );

    /**
    * @brief Scores Target-Decoy pairs using parallel processing.
    *
    * This method scores Target-Decoy pairs based on specified parameters such as the number of top MS2 ions,
    * scan time range, MS calibration information, map of MzTargetKeys to candidate pointers,
    * and a vector to store the resulting candidate scores. The scoring is performed in parallel,
    * leveraging available system processors for efficiency.
    *
    * @param topNMS2Ions Number of top MS2 ions to consider for scoring.
    * @param scanTimeMinMax Pair representing the minimum and maximum scan times for scoring.
    * @param msCalibratomatic MS calibration information.
    * @param mzTargetKeyVsTargetDecoyCandidatePointers Map of MzTargetKeys to vectors of Target-Decoy candidate pointers.
    * @param candidateScoresVec Vector to store the resulting candidate scores.
    * @return Error code indicating success or failure during the scoring process.
    */
    Err scoreTargetDecoyPairs(
            int topNMS2Ions,
            const QPair<ScanTime, ScanTime> &scanTimeMinMax,
            const MsCalibratomatic &msCalibratomatic,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<CandidateScores> *candidateScoresVec
            );

    /**
    * @brief Checks if the TargetDecoyCandidatePairScoretron is initialized.
    *
    * This method checks whether the necessary components for scoring, including Pythia parameters,
    * DIA target frames, and MS1 frame, are initialized.
    *
    * @return True if initialized, false otherwise.
    */
    bool isInit();

    /**
    * @brief Sets the Pythia parameters for TargetDecoyCandidatePairScoretron.
    *
    * This method sets the Pythia parameters and checks if they are valid.
    *
    * @param pythiaParameters The Pythia parameters to set.
    * @return An error code indicating success or failure.
    */
    Err setPythiaParameters(const PythiaParameters &pythiaParameters);

private:

    Err buildParallelInput(
            int topNMS2Ions,
            const QPair<double, double> &scanTimeMinMax,
            const MsCalibratomatic &msCalibratomatic,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<TargetDecoyPairParallelInput> *input
            );

private:

    PythiaParameters m_pythiaParameters;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    TurboXIC *m_turboXICMS1;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *m_diaTargetFrames;
    QMap<ScanNumber, ScanPoints> m_ms1Frame;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
