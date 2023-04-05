//
// Created by anichols on 4/4/23.
//

#ifndef PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H
#define PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H

#include "CSVReader.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "WorkFlowsLib_Exports.h"

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
            QString *outputFilePath
            );

private:

    Err digestFastaEntries(
            const QMap<ProteinId, FastaEntry> &fastaEntries,
            QVector<QSharedPointer<CSVReaderBase>> *peptideSequences
            );

private:

    PythiaParameters m_params;


};


#endif //PYTHIADIACPP_FASTAFILETOPEPTIDESLISTWORKFLOW_H
