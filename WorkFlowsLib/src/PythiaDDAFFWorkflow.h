//
// Created by andrewnichols on 5/22/25.
//

#ifndef PYTHIADDAFFWORKFLOW_H
#define PYTHIADDAFFWORKFLOW_H

#include "WorkFlowsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"

using namespace Error;

class WORKFLOWSLIB_EXPORTS PythiaDDAFFWorkflow {

public:

    PythiaDDAFFWorkflow();
    ~PythiaDDAFFWorkflow();

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
            const QString &fastaUri,
            const QString &outputFolderPath
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

    PythiaParameters m_pythiaParameters;

};



#endif //PYTHIADDAFFWORKFLOW_H
