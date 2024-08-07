//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;


class MsCalibratomatic;
class TargetDecoyPairParallelInput;
class TurboXIC;

class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairScoretron2 {

public:

    TargetDecoyCandidatePairScoretron2();
    ~TargetDecoyCandidatePairScoretron2();

    /**
    * @brief Initializes the TargetDecoyCandidatePairScoretron2.
    *
    * This method initializes the TargetDecoyCandidatePairScoretron2 with the provided parameters,
    * including Pythia parameters, MS reader accessor
    *
    * @param pythiaParameters The Pythia parameters used for scoring.
    * @param msReaderPointerAcc Pointer to the MS reader accessor for accessing MS data.
    * @return Error code indicating success or failure during initialization.
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            MsReaderPointerAcc *msReaderPointerAcc
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
            const MsCalibratomatic &msCalibratomatic,
            float minPeakCount,
            int threadCount,
            const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<CandidateScores> *candidateScoresVec
            ) const;

    /**
    * @brief Checks if the TargetDecoyCandidatePairScoretron2 is initialized.
    *
    * This method checks whether the necessary components for scoring, including Pythia parameters,
    * DIA target frames, and MS1 frame, are initialized.
    *
    * @return True if initialized, false otherwise.
    */
    bool isInit() const;

    /**
    * @brief Sets the Pythia parameters for TargetDecoyCandidatePairScoretron2.
    *
    * This method sets the Pythia parameters and checks if they are valid.
    *
    * @param pythiaParameters The Pythia parameters to set.
    * @return An error code indicating success or failure.
    */
    Err setPythiaParameters(const PythiaParameters &pythiaParameters);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>>* diaTargetFrames();

    QMap<ScanNumber, ScanPoints>* ms1ScanNumberVsScanPoints();

    Err reloadTurboXICMS1();

private:

    Err buildParallelInput(
            int topNMS2Ions,
            const QPair<double, double> &scanTimeMinMax,
            const MsCalibratomatic &msCalibratomatic,
            float minPeakCount,
            const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
            const QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<TargetDecoyPairParallelInput> *input
            ) const;

    Err buildMzTargetKeyVsMsFrames();

    Err buildAveragineTable();


private:

    PythiaParameters m_pythiaParameters;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    TurboXIC *m_turboXICMS1;
    MsFrame *m_msFrameMS1;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> m_diaTargetFrames;
    QMap<ScanNumber, ScanPoints> m_ms1ScanNumberVsScanPoints;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;
    QPair<ScanTime, ScanTime> m_scanTimeMinMax;
    QVector<MsScanInfo> m_uniqueTandemMsScanInfos;
    QMap<MzTargetKey, MsFrame*> m_mzTargetKeyVsMsFrame;
    QMap<NominalMzMass, QVector<float>> m_averagineTable;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
