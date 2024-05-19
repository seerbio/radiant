//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron2.h"

using namespace Error;

class CandidateScores;
class MsScanInfo;
class MsReaderPointerAcc;
class TargetDecoyCandidatePair;

using CandidateScoresTarget = CandidateScores;
using CandidateScoresDecoy = CandidateScores;

struct KarnnNNTarget {
    QString seq;
    float nnScore = 0.0;
    bool isDecoy = false;
    QVector<float> scoreVec;
    int index = -1;
};

class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflow {

public:

    friend class PythiaDIAFFWorkflowTests;

    PythiaDIAFFWorkflow();
    ~PythiaDIAFFWorkflow() = default;

    /**
    * @brief Initializes the PythiaDIAFFWorkflow with necessary parameters
    *
    * @param pythiaParameters Parameters specific to the Pythia engine
    * @param fragLibUri File path for the fragmentation library
    * @param fastaUri File path for the fasta data
    *
    * This function is responsible for setting up the Pythia Data-Independent Acquisition (DIA) workflow.
    * It performs tasks such as validating the provided parameters and files, and calculating bounding mass ranges
    * based on the provided peptide minimum and maximum lengths.
    * After completing the validation tasks, it reads rows from the fragmentation library and initializes the target-decoy candidate pair manager.
    *
    * @return Err Error status of the function execution, which includes any errors occurred during the initialization
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &fastaUri
            );

    /**
    * @brief Executes data analysis and interpretation workflow for MS/MS data
    *
    * @param msDataFilePath File path for the input Mass Spectrometry (MS) data
    *
    * This function is responsible for processing the given Mass Spectrometry data file. It performs a series of steps including:
    * 1. Validating the availability of the input file
    * 2. Gathering relevant data points from the file
    * 3. Building a calibration model
    * 4. Performing main analysis and creating a score for candidates based on the analysis
    * 5. Filtering the candidates based on a certain False Discovery Rate (FDR)
    * 6. Removing interfering candidates
    * 7. Applying a neural network classifier to further refine the candidates
    * 8. Updating Protein Group Annotations
    *
    * The function uses a number of sub-functions and helper classes to achieve the steps. The process includes multiple phases of optimizations as well.
    *
    * @return Err Error status of the function execution, which includes any errors occurred during the file processing
    */
    Err processFile(const QString &msDataFilePath);


private:

    Err buildCalibration(
            MsReaderPointerAcc *msReaderPointerAcc,
            QVector<CandidateScores*> *candidateScoresForTrainings
            );

    Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
            const QVector<MsScanInfo> &msScanInfos,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
            );

    Err recalibrateMzVals(
        const MSLevelEnum &msLevel,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        );

    Err recalibrateMs1Points(const QVector<CandidateScores*> &candidateScoresVecBatchPntrsResized);

    Err optimizeParameters(
            const QVector<CandidateScores*> &candidateScoresTrainings,
            MsReaderPointerAcc *msReaderPointerAcc
            );

    Err mainAnalysis(
            MsReaderPointerAcc *msReaderPointerAcc,
            int *targetCountBelowFDRThresholdOnePercent
            );

    Err buildCandidateScoresPtrs(QVector<CandidateScores*> *candidateScoresPntrs);

    static Err buildCandidateScoresPtrs(
            QVector<CandidateScores> &candidateScores,
            QVector<CandidateScores*> *candidateScoresPntrs
            );

    Err applyNeuralNetClassifier(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            int seed,
            QVector<CandidateScores*> *candidateScoreClassifier
            );

    Err updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<CandidateScores*> *candidateScores
        );

    Err buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(const QVector<CandidateScores*> &candidateScores);

    Err honeIRTAndMassCalibration(
        QVector<CandidateScores*> *candidateScoresVecScoredPntrs,
        int topNCandidates
    );


private:

    TargetDecoyCandidatePairManager m_targetDecoyCandidatePairManager;
    TargetDecoyCandidatePairScoretron2 m_targetDecoyCandidatePairScoretron;
    MsCalibratomatic m_msCalibratomatic;

    QVector<TargetDecoyCandidatePair*> m_targetDecoyPairPntrs;
    QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsTopScores;
    QHash<TargetDecoyCandidatePair*, bool> m_entered;

    QVector<CandidateScores> m_candidateScores;
    QHash<PeptideSequenceWithModsChargeAndTargetKey , QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> m_peptideKeyVsTargetDecoyCandidateScoresPntrs;

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;

    int m_minTopNMs2Ions;
    int m_minTrainingCountTranche;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
