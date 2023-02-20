//
// Created by anichols on 2/17/23.
//

#include "LibraryBuilderWorkFlow.h"

#include "FragLibraryTron.h"
#include "MathUtils.h"
#include "Molecule.h"
#include "PeptidesLibraryTron.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>

namespace {

    QVector<double> vectorizePeptideSequenceAminoAcidMasses(
            const QString &peptideSequence,
            const AminoAcids &aminoAcids
            ) {

        QVector<double> aaMassVals(peptideSequence.size());

        for (int i = 0; i < peptideSequence.size(); i++) {

            const QChar &aa = peptideSequence.at(i);
            const double aaMass = aminoAcids.aminoAcid(aa).monoisotopicMass();
            aaMassVals[i] = aaMass;
        }

        return aaMassVals;
    }

    void addModificationsToVec(
            const QHash<ResidueIndex, ModificationMass> &modifications,
            QVector<double> *vec
            ) {

        for (auto it = modifications.begin(); it != modifications.end(); it++) {

            const ResidueIndex residueIndex = it.key();
            const ModificationMass modificationMass = it.value();

            (*vec)[residueIndex] += modificationMass;
        }
    }

    QVector<double> calcMzPeptideFragmentSeries(
            const Peptide &peptide,
            const AminoAcids &aminoAcids,
            bool isYSeries
            ) {

        QVector<double> aaMassVals = vectorizePeptideSequenceAminoAcidMasses(peptide.sequence, aminoAcids);

        addModificationsToVec(
                peptide.modifications,
                &aaMassVals
                );

        aaMassVals[aaMassVals.size() - 1] += CommonMolecules::H2O.monoisotopicMass();

        if (isYSeries) {
            std::reverse(aaMassVals.begin(), aaMassVals.end());
        }

        aaMassVals[0] += PROTON;

        QVector<double> cumSumSeries(aaMassVals.size());
        std::partial_sum(aaMassVals.begin(), aaMassVals.end(), cumSumSeries.begin(), std::plus<>());

        return cumSumSeries;
    }

    void filterMzMasses(
            double mzMin,
            double mzMax,
            QVector<double> *vec
            ) {

        const auto terminatorLogic = [mzMin, mzMax](double mz){
            return !(mzMin <= mz && mz <= mzMax);
        };

        const auto terminator = std::remove_if(vec->begin(), vec->end(), terminatorLogic);
        vec->erase(terminator, vec->end());
    }

    QVector<double> buildFragSeriesBY(
            const Peptide &peptide,
            const AminoAcids &aminoAcids,
            double mzMin,
            double mzMax
            ) {

        QVector<double> bySeries;

        QVector<double> bSeries = calcMzPeptideFragmentSeries(
                peptide,
                aminoAcids,
                false
                );
        bSeries.pop_back();

        QVector<double> ySeries = calcMzPeptideFragmentSeries(
                peptide,
                aminoAcids,
                true
        );
        ySeries.pop_back();

        bySeries.append(bSeries);
        bySeries.append(ySeries);

        const int bySeriesSize = bySeries.size();
        for (int i = 0; i < bySeriesSize; ++i) {
            const double fragIonCharge2 = (bySeries.at(i) + PROTON) / 2;
            bySeries.push_back(fragIonCharge2);
        }

        std::sort(bySeries.begin(), bySeries.end());
        filterMzMasses(
                mzMin,
                mzMax,
                &bySeries
                );

        return bySeries;
    }

    struct ParallelFragInput {
        Peptide peptide;
        PythiaParameters params;
    };

    QVector<ParallelFragInput> buildParallelFragInputs(
            const QVector<Peptide> &peptides,
            const PythiaParameters &params
            ) {

        QVector<ParallelFragInput> inputs;
        for (const Peptide &peptide : peptides) {

            ParallelFragInput pi;
            pi.peptide = peptide;
            pi.params = params;

            inputs.push_back(pi);
        }

        return inputs;
    }

    QPair<Peptide, QVector<double>> parallelFragLogic(const ParallelFragInput &pfi) {
        return {pfi.peptide,
                buildFragSeriesBY(
                pfi.peptide,
                pfi.params.aminoAcids,
                pfi.params.mzMinDataStructure,
                pfi.params.mzMaxDataStructure
                )};
    }

    Err writeFragLibIons(
            const QVector<QPair<Peptide, QVector<double>>> &mzFrags,
            const QString &fragLibIonsFilePath
            ) {

        ERR_INIT

        QVector<FragLibIon> fragLibIons;
        for (const QPair<Peptide, QVector<double>> &qp : mzFrags) {

            const Peptide &peptide = qp.first;
            const QVector<double> &fragIons = qp.second;

            for (double mz : fragIons) {
                fragLibIons.push_back({peptide.id, mz, peptide.mass});
            }
        }

        e = FragLibraryTron::writeFragLibIons(
                fragLibIons,
                fragLibIonsFilePath
                ); ree;

        ERR_RETURN
    }


}//namespace
Err LibraryBuilderWorkFlow::exec(
        const PythiaParameters &pythiaParameters,
        const QString &fastaFilePath,
        bool theoreticalFrag
        ) {

    ERR_INIT

    m_pythiaParameters = pythiaParameters;

    PeptidesLibraryTron peptidesLibraryTron(pythiaParameters);
    e = peptidesLibraryTron.exec(fastaFilePath, 666); ree;

    const QVector<Peptide> peptides = peptidesLibraryTron.peptides();
    e = ErrorUtils::isNotEmpty(peptides);

    QVector<QPair<Peptide, QVector<double>>> mzFrags;

    if (theoreticalFrag) {

        e = buildTheoreticalMzFragsForPeptides(
                peptides,
                &mzFrags
                ); ree;
    }

    else {
        //TODO add neural net frag predictions
    }

    const QString fragLibIonsFilePath = fastaFilePath + S_GLOBAL_SETTINGS.DOT_FRAGLIB;
    e = writeFragLibIons(
            mzFrags,
            fragLibIonsFilePath
            ); ree;

    const QString peptidesLibraryTronFilePath = fastaFilePath + S_GLOBAL_SETTINGS.DOT_PEPLIB;
    e = PeptidesLibraryTron::writePeptidesLib(
            peptides,
            peptidesLibraryTronFilePath
            ); ree;

    ERR_RETURN
}

Err LibraryBuilderWorkFlow::buildTheoreticalMzFragsForPeptides(
        const QVector<Peptide> &peptides,
        QVector<QPair<Peptide, QVector<double>>> *mzFrags
        ) {

    ERR_INIT

#define RUN_PARALLEL_THEO_FRAGMENTATION
#ifdef RUN_PARALLEL_THEO_FRAGMENTATION
    const QVector<ParallelFragInput> inputs = buildParallelFragInputs(
            peptides,
            m_pythiaParameters
    );

    QFuture<QPair<Peptide, QVector<double>>> futures = QtConcurrent::mapped(
            inputs,
            parallelFragLogic
    );
    futures.waitForFinished();

    for (const QPair<Peptide, QVector<double>> &res : futures) {
        mzFrags->push_back(res);
    }
#else
    for (const Peptide &pep: peptides) {

        const QVector<double> bySeries = buildFragSeriesBY(
                pep,
                m_pythiaParameters.aminoAcids,
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure
        );


        mzFrags->push_back({pep, bySeries});
    }
#endif

    ERR_RETURN
}

QVector<double> LibraryBuilderWorkFlow::testPeptideFragmentation(
        const QString &peptideSequence,
        const QHash<ResidueIndex, ModificationMass> &mods
        ) {

    Peptide pep;
    pep.sequence = peptideSequence;
    pep.modifications =mods;

    return buildFragSeriesBY(pep, AminoAcids(), 0.0, 2000.0);
}
