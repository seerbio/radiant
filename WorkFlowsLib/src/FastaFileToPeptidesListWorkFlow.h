//
// Created by anichols on 4/4/23.
//

#ifndef PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H
#define PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H

#include "ParquetReader.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "WorkFlowsLib_Exports.h"
#include "PeptideStringWithMods.h"

#include "Error.h"

using namespace Error;

class FastaEntry;
class PeptideSequence;

class WORKFLOWSLIB_EXPORTS FastaFileToPeptidesListWorkFlow {

public:

    FastaFileToPeptidesListWorkFlow() = default;
    ~FastaFileToPeptidesListWorkFlow() = default;

    Err init(const PythiaParameters &pythiaParameters);

    Err exec(
            const QString &fastaFilePath,
            const QString &targetMzCollisionCSVFilePath,
            const QString &outputFilePath
            );

    static Err digestFastaEntries(
            const PythiaParameters &pythiaParameters,
            const QMap<ProteinId, FastaEntry> &fastaEntries,
            QVector<PeptideSequence> *peptideSequences,
            QMap<PeptideStringWithMods, QVector<FastaEntry>> *peptideStringWithModsVsFastaEntries
    );

private:

    Err digestFastaEntries(
            const QMap<ProteinId, FastaEntry> &fastaEntries,
            QVector<PeptideSequence> *peptideSequences
    );

    Err addDecoys(
            int seed,
            QVector<PeptideSequence> *peptideSequences
    );

private:

    PythiaParameters m_params;


};


#endif //PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H
