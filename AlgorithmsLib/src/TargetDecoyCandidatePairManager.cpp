//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePairManager.h"

#include "AminoAcids.h"
#include "FragLibReader.h"
#include "MolecularFormula.h"
#include "Molecule.h"

#include <QElapsedTimer>

namespace {

    Err buildMzIons(
            const FragLibReaderRow &row,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(row.mzVals); ree;
        e = ErrorUtils::isEqual(row.mzVals.size(), row.intensityVals.size()); ree;

        const QStringList ionLabels = row.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

        e = ErrorUtils::isEqual(row.mzVals.size(), ionLabels.size()); ree;

        for (int i = 0; i < row.mzVals.size(); i++) {
            MS2Ion ion;
            ion.mz = row.mzVals.at(i);
            ion.intensity = row.intensityVals.at(i);
            ion.ionLabel = ionLabels.at(i);
            ion.iRT = static_cast<float>(row.iRT);
            ion.charge = ion.ionLabel.contains("^2") ? 2 : 1;
            ms2Ions->push_back(ion);
        }

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairManager::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibFileUri
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibFileUri); ree;

    m_pythiaParameters = pythiaParameters;

    const double massMin
        = pythiaParameters.peptideLengthMin * Molecule(MolecularFormulas::alanineFormula).monoisotopicMass();

    const double massMax
            = pythiaParameters.peptideLengthMax * Molecule(MolecularFormulas::tryptophanFormula).monoisotopicMass();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            fragLibFileUri,
            massMin,
            massMax,
            &fragLibReaderRows
    );

    QElapsedTimer et;
    et.start();

    e = buildTargetDecoyCandidatePairs(fragLibReaderRows); ree;
    e = buildIndexVsTargetDecoyCandidatePairPtrs(); ree;

    qDebug() << m_targetDecoyCandidatePairs.size() << "Candidates loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    Err mutateCandidatePeptideTarget(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonTarget,
            QVector<MS2Ion> *ms2IonDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;
        e = ErrorUtils::isNotEmpty(ms2IonTarget); ree;

        ms2IonDecoys->clear();
        ms2IonDecoys->reserve(ms2IonTarget.size());

        const QMap<QChar, double> diannMutateAminoAcidTo = AminoAcids::diannMutateAminoAcidTo();

        const QString &seq = peptideStringWithMods;

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = seq.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }
                }
            }

            else if (ionLableInfo.second.contains('y')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }

                    if (ionLableInfo.first >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }
                }
            }

            else {
                qDebug() << "Non b/y/a ion" << ionLableInfo;
            }

            ms2IonDecoys->push_back(ms2IonDecoy);
        }

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairManager::buildTargetDecoyCandidatePairs(
        const QVector<FragLibReaderRow> &fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibReaderRows); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;

    m_targetDecoyCandidatePairs.reserve(fragLibReaderRows.size()); ree;

    for (const FragLibReaderRow &flrr : fragLibReaderRows) {

        PeptideStringWithMods peptideStringWithMods;
        Charge charge;
        e = peptideStringWithModsFromPeptideSequenceChargeKey(
                flrr.peptideSequenceChargeKey,
                &peptideStringWithMods,
                &charge
        ); ree;

        //TODO once incorporating PTMs in the sequence, make algo to remove them
        // and then calculate the peptide length MODS
        const int peptideLength = peptideStringWithMods.size();
        if (peptideLength < m_pythiaParameters.peptideLengthMin
            || peptideLength > m_pythiaParameters.peptideLengthMax
            || flrr.isDecoy) {
            continue;
        }

        QVector<MS2Ion> ms2IonsTarget;
        e = buildMS2Ions(flrr, &ms2IonsTarget); ree;

        QVector<MS2Ion> ms2IonsDecoy;
        e= mutateCandidatePeptideTarget(
                peptideStringWithMods,
                ms2IonsTarget,
                &ms2IonsDecoy
                ); ree;

        TargetDecoyCandidatePair targetDecoyCandidatePair(
                peptideStringWithMods,
                ms2IonsTarget,
                ms2IonsDecoy,
                flrr.charge,
                flrr.mass,
                flrr.iRT,
                flrr.mzVals.size(),
                m_targetDecoyCandidatePairs.size()
                );

        m_targetDecoyCandidatePairs.push_back(targetDecoyCandidatePair);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::buildMS2Ions(
        const FragLibReaderRow &flrr,
        QVector<MS2Ion> *ms2Ions
) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(flrr.mzVals); ree;
    e = ErrorUtils::isEqual(flrr.mzVals.size(), flrr.intensityVals.size());ree;
    e = ErrorUtils::isNotEmpty(flrr.ionLabels); ree

    const QStringList ionLabelsSplit = flrr.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);
    e = ErrorUtils::isEqual(flrr.mzVals.size(), ionLabelsSplit.size());ree;

    QVector<MS2Ion> ms2IonsBuilder;
    ms2IonsBuilder.reserve(flrr.mzVals.size());

    for (int i = 0; i < flrr.mzVals.size(); i++) {
        MS2Ion ms2Ion;
        ms2Ion.mz = flrr.mzVals.at(i);
        ms2Ion.intensity = flrr.intensityVals.at(i);
        ms2Ion.ionLabel = ionLabelsSplit.at(i);
        ms2Ion.iRT = static_cast<IRT>(flrr.iRT);
        ms2Ion.charge = flrr.charge;

        ms2IonsBuilder.push_back(ms2Ion);
    }

    MS2Ion::filterMS2IonsByMz(
            m_pythiaParameters.mzMinDataStructure,
            m_pythiaParameters.mzMaxDataStructure,
            &ms2IonsBuilder
            );

    MS2Ion::sortMS2IonsIntensityDesc(&ms2IonsBuilder);

    const int ms2IonsSize = std::min(m_pythiaParameters.topNMs2Ions, ms2IonsBuilder.size());
    ms2Ions->reserve(ms2IonsSize);
    ms2IonsBuilder.resize(ms2IonsSize);

    *ms2Ions = ms2IonsBuilder;

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::buildIndexVsTargetDecoyCandidatePairPtrs() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;
    for (TargetDecoyCandidatePair &tdcp : m_targetDecoyCandidatePairs) {
        m_indexVsTargetDecoyCandidatePairPtrs.insert(tdcp.targetDecoyCandidatePairIndex(), &tdcp);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        PeptideStringWithMods *peptideStringWithMods,
        Charge *charge
){

    ERR_INIT

    const int expectedSplitSize = 2;

    const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
            S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
            QString::SkipEmptyParts
    );

    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize); ree;

    *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
    ); ree

    ERR_RETURN
}

