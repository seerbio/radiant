//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
#define PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


using namespace Error;


class Peptide;

class WORKFLOWSLIB_EXPORTS LibraryBuilderWorkFlow {

public:

    friend class LibraryBuilderWorkFlowTests;

    LibraryBuilderWorkFlow() = default;
    ~LibraryBuilderWorkFlow() = default;

    Err exec(
            const PythiaParameters &pythiaParameters,
            const QString &fastaFilePath,
            bool theoreticalFrags
            );

private:

    Err buildTheoreticalMzFragsForPeptides(
            const QVector<Peptide> &peptides,
            QVector<QPair<Peptide, QVector<double>>> *mzFrags
            );

    static QVector<double> testPeptideFragmentation(
            const QString &peptideSequence,
            const QHash<ResidueIndex, ModificationMass> &mods
            );

private:

    PythiaParameters m_pythiaParameters;


};


#endif //PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
