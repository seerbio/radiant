//
// Created by anichols on 2/16/23.
//

#include "PeptidesLibraryTron.h"

#include "BiophysicalCalcs.h"
#include "FastaReader.h"
#include "Molecule.h"
#include "MolecularFormula.h"
#include "StringUtils.h"

#include <QElapsedTimer>
#include <QList>
#include <QtConcurrent/QtConcurrent>

#include <iostream>
#include <random>


PeptidesLibraryTron::PeptidesLibraryTron(const PythiaParameters &pythiaParameters)
: m_pythiaParameters(pythiaParameters) {}

namespace {

    int countDecoys(const QVector<Peptide> &peptides) {

        int counter = 0;
        for (const Peptide &p : peptides) {
            if (p.isDecoy) {
                counter++;
            }
        }

        return counter;
    }

}//namespace
Err PeptidesLibraryTron::exec(
        const QString &fastaFileUri,
        int seed
        ){

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = processFasta(fastaFileUri); ree;
    e = removeDuplicatesFromPeptide(); ree;
    e = addDecoys(seed); ree;
    e = buildPeptides(); ree;
    e = addVariableModificationsToPeptides(); ree;
    e = addTerminalModificationsToPeptides(); ree;
    e = addPeptideIdToPeptides(); ree;
    e = addMassToPeptides(); ree;

    const int decoyCount = countDecoys(m_peptides);
    const int targetCount = m_peptides.size() - decoyCount;

    qDebug() << "Finished in" << et.elapsed() << "mSec";
    qDebug() << "Total Peptide Count" << m_peptides.size();
    qDebug() << "targets:" << targetCount << "decoys:" << decoyCount;

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

    // Checks to see if the modification is specific to
    // the N or C terminal of the Protein as a whole, or
    // individual peptide.
    bool isTerminalModification(const Modification &mod) {
        return !mod.positionalLocation.isEmpty();
    }

    QPair<Err,QVector<Peptide>> iterateModification(
            const Peptide &peptide,
            const PythiaParameters &pythiaParameters
            ) {

        ERR_INIT

        QVector<Peptide> modPeptides;

        for (const Modification &mod : pythiaParameters.modifications) {


            if (mod.type == ModificationType::FIXED || isTerminalModification((mod))) {
                continue;
            }

            MolecularFormula mf;
            e = parseMolecularFormulaString(mod.formula, &mf);
            if (e != eNoError) {
                return {e, {}};
            }

            Molecule mol(mf);

            const QVector<int> indexesOfModResidue = StringUtils::findIndexesOfCharacterInString(
                    peptide.sequence,
                    mod.residue
                    );

            for (int residueIndex : indexesOfModResidue) {

                if (peptide.modifications.keys().contains(residueIndex)) {
                    continue;
                }

                Peptide newModPeptide = peptide;
                newModPeptide.modifications.insert(residueIndex, mol.monoisotopicMass());

                modPeptides.push_back(newModPeptide);
            }

        }

        return {e, modPeptides};
    }

    QPair<Err,QVector<Peptide>> addVariableModToPeptide(
            const Peptide &peptide,
            const PythiaParameters &pythiaParameters
            ) {

        QPair<Err,QVector<Peptide>> peptidesModified = iterateModification(
                peptide,
                pythiaParameters
                );
        if (peptidesModified.first != eNoError) {
            return {peptidesModified.first, {}};
        }

        QMap<QString, bool> entered;

        for (int maxModIndexCounter = 1;
                maxModIndexCounter < pythiaParameters.maxModificationsPeptide;
                maxModIndexCounter++) {

            const int peptidesMofifiedSize = peptidesModified.second.size();

            for (int i = 0; i < peptidesMofifiedSize; i++) {

                const Peptide &pepMod = peptidesModified.second[i];

                if (pepMod.modifications.size() > pythiaParameters.maxModificationsPeptide) {
                    continue;
                }

                const QPair<Err,QVector<Peptide>> newModPeptides = iterateModification(pepMod, pythiaParameters);
                if (peptidesModified.first != eNoError) {
                    return {peptidesModified.first, {}};
                }

                for (const Peptide &np : newModPeptides.second) {

                    const QString enteredKey = np.peptideStringWithMods();
                    if (entered.value(enteredKey)) {
                        continue;
                    }
                    peptidesModified.second.push_back(np);
                    entered.insert(enteredKey, true);
                }
            }
        }

        return peptidesModified;
    }

    struct InputParallelVarMods {
        Peptide peptide;
        PythiaParameters params;
    };

    QVector<InputParallelVarMods> buildInputeParallelVarMods(
            const QVector<Peptide> &peptides,
            const PythiaParameters &params
            ) {

        QVector<InputParallelVarMods> inputs;

        for (const Peptide &peptide : peptides) {
            InputParallelVarMods inputParallelVarMod;
            inputParallelVarMod.peptide = peptide;
            inputParallelVarMod.params = params;

            inputs.push_back(inputParallelVarMod);
        }

        return inputs;
    }

    QPair<Err,QVector<Peptide>> addVariableModParallelLogic(const InputParallelVarMods &input) {
        return addVariableModToPeptide(input.peptide, input.params);
    }

}//namespace
Err PeptidesLibraryTron::addVariableModificationsToPeptides() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptides); ree;

#define RUN_PARALLEL_ADD_VAR_MODS
#ifdef RUN_PARALLEL_ADD_VAR_MODS

    QVector<InputParallelVarMods> inputs = buildInputeParallelVarMods(
            m_peptides,
            m_pythiaParameters
            );

    QFuture<QPair<Err,QVector<Peptide>>> futures = QtConcurrent::mapped(
            inputs,
            addVariableModParallelLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<Peptide>> & result : futures) {

        if (result.first != eNoError) {
            rrr(eError);
        }

        m_peptides.append(result.second);
    }

#else
    QVector<Peptide> moddedPeptides;

    for (const Peptide &pep : m_peptides) {

        QPair<Err,QVector<Peptide>> modPeptides = addVariableModToPeptide(pep, m_pythiaParameters);
        if (modPeptides.first != eNoError) {
            rrr(eError);
        }
        moddedPeptides.append(modPeptides.second);
    }

    m_peptides.append(moddedPeptides);
#endif

    ERR_RETURN
}

Err PeptidesLibraryTron::addTerminalModificationsToPeptides() {

    ERR_INIT

    for (const Modification &mod : m_pythiaParameters.modifications) {

        if (mod.type == ModificationType::FIXED || !isTerminalModification((mod))) {
            continue;
        }

        MolecularFormula mf;
        e = parseMolecularFormulaString(mod.formula, &mf); ree;
        Molecule mol(mf);

        const int currentPeptidesSize = m_peptides.size();
        for (int i = 0; i < currentPeptidesSize; i++) {

            Peptide newPeptide = m_peptides.at(i);

            if (mod.positionalLocation == Modification::nTermProtein() && newPeptide.previousResidue == '*') {

                newPeptide.modifications.insert(0, mol.monoisotopicMass());

            } else if (mod.positionalLocation == Modification::cTermProtein() && newPeptide.postResidue == '*') {

                newPeptide.modifications.insert(newPeptide.sequence.size() - 1, mol.monoisotopicMass());
            }

            //TODO add logic for peptide terminal modifications.

            else {
                continue;
            }

            m_peptides.push_back(newPeptide);
        }
    }

    ERR_RETURN
}

Err PeptidesLibraryTron::addPeptideIdToPeptides() {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_peptides); ree;
    for (int i = 0; i < m_peptides.size(); i++) {
        Peptide &p = m_peptides[i];
        p.id = i;
    }

    e = ErrorUtils::isTrue(m_peptides.last().id == m_peptides.size() - 1); ree;

    ERR_RETURN
}

Err PeptidesLibraryTron::addMassToPeptides() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_peptides); ree;

    for (int i = 0; i < m_peptides.size(); i++) {
        Peptide &pep = m_peptides[i];
        pep.mass = BiophysicalCalcs::calculatePeptideMass(
                pep.sequence,
                m_pythiaParameters.aminoAcids,
                pep.modifications
                );
    }

    e = ErrorUtils::isTrue(m_peptides.first().mass > 0);
    e = ErrorUtils::isTrue(m_peptides.last().mass > 0);

    ERR_RETURN
}

QVector<Peptide> PeptidesLibraryTron::peptides() const {
    return m_peptides;
}
