//
// Created by anichols on 2/16/23.
//

#include "PeptidesLibraryTron.h"

#include "FastaReader.h"
#include "Molecule.h"
#include "MolecularFormula.h"
#include "StringUtils.h"

#include <QElapsedTimer>
#include <QList>
#include <QtAlgorithms>
#include <QtConcurrent/QtConcurrent>

#include <iostream>
#include <random>


PeptidesLibraryTron::PeptidesLibraryTron(const PythiaParameters &pythiaParameters)
: m_pythiaParameters(pythiaParameters)
{}

Err PeptidesLibraryTron::exec(
        const QString &fastaFileUri,
        int seed
        ){

    ERR_INIT

    e = processFasta(fastaFileUri); ree;
    e = removeDuplicatesFromPeptide(); ree;
    e = addDecoys(seed); ree;
    e = buildPeptides(); ree;
    e = addVariableModificationsToPeptides(); ree;

    ERR_RETURN
}

namespace {

    struct ParallelProcessFastaInput {
        PythiaParameters pythiaParameters;
        FastaEntry fastaEntry;
    };

    QVector<ParallelProcessFastaInput> buildParallelProcessFastaInput(
            const PythiaParameters &pythiaParameterReader,
            const QMap<ProteinId, FastaEntry> &fastaEntries
            ) {

        QVector<ParallelProcessFastaInput> output;

        for (const FastaEntry &fe : fastaEntries) {

            ParallelProcessFastaInput pi;
            pi.pythiaParameters = pythiaParameterReader;
            pi.fastaEntry = fe;
            output.push_back(pi);
        }

        return output;
    }

    QVector<PeptideSequence> parallelProcessFastaLogic(const ParallelProcessFastaInput &input) {

        ProteinDigestomatic proteinDigestomatic(input.pythiaParameters);

        QVector<PeptideSequence> peptideSequences;

        proteinDigestomatic.digestProtein(
                input.fastaEntry.fastaSequence,
                &peptideSequences
                );

        return peptideSequences;
    }

}//namespace
Err PeptidesLibraryTron::processFasta(const QString &fastaFileUri) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fastaFileUri); ree;

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFileUri); ree;
    const QMap<ProteinId, FastaEntry> fastaEntries = fastaReader.fastaEntries();

    e = ErrorUtils::isNotEmpty(fastaEntries); ree;

    QVector<ParallelProcessFastaInput> parallelProcessingInputs = buildParallelProcessFastaInput(
            m_pythiaParameters,
            fastaEntries
            );

    e = ErrorUtils::isNotEmpty(parallelProcessingInputs); ree;

    QFuture<QVector<PeptideSequence>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProcessFastaLogic
            );
    futures.waitForFinished();

    for (const QVector<PeptideSequence> & pepSeqs : futures) {
        m_peptideSequences.append(pepSeqs);
    }

    e = ErrorUtils::isNotEmpty(m_peptideSequences); ree;

    ERR_RETURN
}

namespace {

    void replaceLeucinesWithX(PeptideSequence *pepSeq) {
        pepSeq->sequence.replace("L", "X");
        pepSeq->sequence.replace("I", "X");
    }

}//namespace
Err PeptidesLibraryTron::removeDuplicatesFromPeptide() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptideSequences); ree;

    QMap<QString, bool> alreadyInserted;
    QVector<PeptideSequence> uniqueSeqs;

    for (const PeptideSequence &peptideSequence : m_peptideSequences) {

        PeptideSequence newPepSeq = peptideSequence;
        replaceLeucinesWithX(&newPepSeq);

        if (alreadyInserted.value(newPepSeq.sequence)) {
            continue;
        }

        uniqueSeqs.push_back(newPepSeq);
        alreadyInserted.insert(newPepSeq.sequence, true);
    }

    m_peptideSequences.clear();
    m_peptideSequences = uniqueSeqs;

    ERR_RETURN
}

namespace {

    const QString FAILED_SHUFFLE = QStringLiteral("Failed Shuffle");

    QString shufflePeptide(const QString& peptideSeq) {

        const int peptideLen = peptideSeq.size();
        if (peptideLen < 3) {
            return peptideSeq;
        }

        QStringList strList = peptideSeq.split("");

        std::shuffle(
                strList.begin() + 2,
                strList.begin() + peptideLen - 1,
                std::mt19937(std::random_device()())
                );

        return strList.join("");
    }

    PeptideSequence parallelAddDecoysLogic(const PeptideSequence &ps) {

        PeptideSequence newPepSeq = ps;
        newPepSeq.isDecoy = true;

        int retries = 0;
        const int maxRetryCount = 10;
        while (newPepSeq.sequence == ps.sequence) {

            newPepSeq.sequence = shufflePeptide(ps.sequence);

            if (retries++ > maxRetryCount) {
                newPepSeq.sequence = FAILED_SHUFFLE;
                return newPepSeq;
            }
        }

        return newPepSeq;
    }

}//namespace
Err PeptidesLibraryTron::addDecoys(int seed) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptideSequences); ree;
    const int peptideSequencesSize = m_peptideSequences.size();

    std::srand(666); //TODO fix seeding reproducibility of the shuffle.

#define RUN_PARALLEL_DECOYS
#ifdef RUN_PARALLEL_DECOYS

    QFuture<PeptideSequence> futures = QtConcurrent::mapped(
            m_peptideSequences,
            parallelAddDecoysLogic
    );
    futures.waitForFinished();

    for (const PeptideSequence &pepSeq : futures) {

        if (pepSeq.sequence == FAILED_SHUFFLE) {
            continue;
        }

        m_peptideSequences.append(pepSeq);
    }

    e = ErrorUtils::isTrue(peptideSequencesSize != m_peptideSequences.size()); ree;

#else
    int counter = 0;
    for (const PeptideSequence &ps : m_peptideSequences) {
        const PeptideSequence newPepSeq = parallelAddDecoysLogic(ps);

        if (newPepSeq.sequence == FAILED_SHUFFLE) {
            continue;
        }

        m_peptideSequences.append(newPepSeq);
    }
#endif

    ERR_RETURN
}

Err PeptidesLibraryTron::buildPeptides() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptideSequences); ree;

    QMap<QString, bool> peptideHasBeenEntered;

    for (const PeptideSequence &ps : m_peptideSequences) {

        if (peptideHasBeenEntered.value(ps.sequence)) {
            continue;
        }

        Peptide pep;
        pep.sequence = ps.sequence;
        pep.isDecoy = ps.isDecoy;
        pep.previousResidue = ps.previousResidue;
        pep.postResidue = ps.postResidue;

        m_peptides.push_back(pep);
        peptideHasBeenEntered.insert(ps.sequence, true);
    }

    ERR_RETURN
}

namespace {

    QVector<Peptide> iterateModification(
            const Peptide &peptide,
            const PythiaParameters &pythiaParameters
            ) {

        QVector<Peptide> modPeptides;

        for (const Modification &mod : pythiaParameters.modifications) {


            if (mod.type == ModificationType::FIXED) {
                continue;
            }

            MolecularFormula mf;
            parseMolecularFormulaString(mod.formula, &mf);
            Molecule mol(mf);

            const QVector<int> indexesOfModResidue = StringUtils::findIndexesOfCharacterInString(
                    peptide.sequence,
                    mod.residue
                    );

            for (int residueIndex : indexesOfModResidue) {

                if (peptide.modifications.contains(residueIndex)) {
                    continue;
                }

                Peptide newModPeptide = peptide;
                newModPeptide.modifications.insert(residueIndex, mol.monoisotopicMass());

                modPeptides.push_back(newModPeptide);
            }

        }

        return modPeptides;
    }

    QVector<Peptide> addVariableModToPeptide(
            const Peptide &peptide,
            const PythiaParameters &pythiaParameters
            ) {


        QVector<Peptide> peptidesModified = iterateModification(
                peptide,
                pythiaParameters
                );

        QMap<QString, bool> entered;

        for (int maxModIndexCounter = 1;
                maxModIndexCounter < pythiaParameters.maxModificationsPeptide;
                maxModIndexCounter++) {

            for (const Peptide &pepMod : peptidesModified) {

                if (pepMod.modifications.size() > pythiaParameters.maxModificationsPeptide) {
                    continue;
                }

                QVector<Peptide> newModPeptides = iterateModification(pepMod, pythiaParameters);
                for (const Peptide &np : newModPeptides) {
                    const QString enteredKey = np.peptideStringWithMods();
                    if (entered.value(enteredKey)) {
                        continue;
                    }
                    peptidesModified.push_back(np);
                    entered.insert(enteredKey, true);
                }
            }
        }

//        for (auto p : peptidesModified) {
//            qDebug() << p.peptideStringWithMods();
//        }

        return peptidesModified;
    }


}//namespace
Err PeptidesLibraryTron::addVariableModificationsToPeptides() {

    ERR_INIT

//    e = ErrorUtils::isNotEmpty(m_peptides); ree;

    Peptide pep;
    pep.sequence = "QMASSMSSNSLKSNDMSK";
    pep.previousResidue = '*';

    addVariableModToPeptide(pep, m_pythiaParameters);


    ERR_RETURN
}
