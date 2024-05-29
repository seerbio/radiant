//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PeptideStringWithMods.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"


using namespace Error;

class FragLibReaderRow;


class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairManager {

public:

    friend class TargetDecoyCandidatePairManagerTests;

    TargetDecoyCandidatePairManager();
    ~TargetDecoyCandidatePairManager();

    /**
    * @brief Initializes the TargetDecoyCandidatePairManager.
    *
    * This method performs the following steps:
    * 1. Checks if the PythiaParameters are valid.
    * 2. Ensures that fragLibReaderRows is not empty.
    * 3. Sets the PythiaParameters for the manager.
    * 4. Removes non-sequences with non-canonical amino acids from fragLibReaderRows.
    * 5. Builds target-decoy candidate pairs based on fragLibReaderRows.
    * 6. Initializes the manager's private implementation with the generated candidate pairs.
    *
    * @param pythiaParameters The PythiaParameters to use for initialization.
    * @param fragLibReaderRows A pointer to the QVector<FragLibReaderRow> containing fragment library reader rows.
    * @return An error code indicating the success or failure of the initialization process.
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            QVector<FragLibReaderRow> *fragLibReaderRows
            );

    /**
    * @brief Checks whether the TargetDecoyCandidatePairManager has been successfully initialized.
    *
    * This method returns true if the list of target-decoy candidate pairs is not empty, indicating a successful initialization.
    *
    * @return True if the manager is initialized; false otherwise.
    */
    bool isInit();

    /**
    * @brief Gets the total count of target-decoy candidate pairs managed by the TargetDecoyCandidatePairManager.
    *
    * @return The total count of target-decoy candidate pairs.
    */
    int targetsCount();

    /**
    * @brief Retrieves pointers to TargetDecoyCandidatePair instances within a specified mz range.
    *
    * This method performs the following steps:
    * 1. Checks if the list of target-decoy candidate pairs is not empty.
    * 2. Ensures that mzMax is greater than or equal to mzMin.
    * 3. Calls the private implementation to get pointers to candidate pairs within the specified mz range.
    * 4. Shuffles the obtained pointers using a random number generator.
    *
    * @param mzMin The minimum mz value for filtering candidate pairs.
    * @param mzMax The maximum mz value for filtering candidate pairs.
    * @param targetDecoyPointers A pointer to a QVector<TargetDecoyCandidatePair*> to store the obtained pointers.
    * @return An error code indicating the success or failure of the operation.
    */
    Err getTargetDecoyCandidatePairPointers(
            double mzMin,
            double mzMax,
            QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
            );

    /**
    * @brief Gets pointers to target-decoy candidate pairs within the specified mz range, with optional random selection.
    *
    * @param mzMin The minimum mz value.
    * @param mzMax The maximum mz value.
    * @param randomSelectionFraction The fraction of randomly selected pairs (use -1 for no random selection).
    * @param targetDecoyPointers Pointer to a QVector to store the selected target-decoy candidate pair pointers.
    * @return Error code indicating success or failure.
    */
    Err getTargetDecoyCandidatePairPointers(
            double mzMin,
            double mzMax,
            double randomSelectionFraction,
            QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
    );

    Err getTargetDecoyCandidatePairPointers(QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairsPntrs);

    /**
    * @brief Extracts peptide string with modifications and charge from a peptide sequence charge key.
    *
    * @param peptideSequenceChargeKey The peptide sequence charge key to extract information from.
    * @param peptideStringWithMods Pointer to store the extracted peptide string with modifications.
    * @param charge Pointer to store the extracted charge.
    * @return Error code indicating success or failure.
    */
    static Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
            );

private:

    Err buildTargetDecoyCandidatePairs(const QVector<FragLibReaderRow> &fragLibReaderRows);
    Err filterDecoySequencesThatAreAlsoTargetSequences();

private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    PythiaParameters m_pythiaParameters;

    Q_DISABLE_COPY(TargetDecoyCandidatePairManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
