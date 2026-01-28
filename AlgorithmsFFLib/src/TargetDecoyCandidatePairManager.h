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

    TargetDecoyCandidatePairManager() = default;
    ~TargetDecoyCandidatePairManager() = default;

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
    * @param fragLibReaderRows A pointer to the QList<FragLibReaderRow> containing fragment library reader rows.
    * @return An error code indicating the success or failure of the initialization process.
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            QList<FragLibReaderRow> *fragLibReaderRows
            );

    /**
    * @brief Checks whether the TargetDecoyCandidatePairManager has been successfully initialized.
    *
    * This method returns true if the list of target-decoy candidate pairs is not empty, indicating a successful initialization.
    *
    * @return True if the manager is initialized; false otherwise.
    */
    bool isInit() const;

    /**
    * @brief Gets the total count of target-decoy candidate pairs managed by the TargetDecoyCandidatePairManager.
    *
    * @return The total count of target-decoy candidate pairs.
    */
    int targetsCount() const;

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
            int *charge
            );

private:

    Err buildTargetDecoyCandidatePairs(QList<FragLibReaderRow> *fragLibReaderRows);
    Err filterDecoySequencesThatAreAlsoTargetSequences();

private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    PythiaParameters m_pythiaParameters;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
