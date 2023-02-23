//
// Created by anichols on 2/16/23.
//

#ifndef PYTHIADIACPP_PEPTIDESLIBRARYTRON_H
#define PYTHIADIACPP_PEPTIDESLIBRARYTRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"

#include <QDataStream>
#include <QStringList>


using namespace Error;

class ALGORITHMSLIB_EXPORTS Peptide {

        public:

        Peptide() = default;
        ~Peptide() = default;

        QString sequence;
        double mass = -1.0;
        PeptideId id = -1;
        bool isDecoy = false;
        QChar previousResidue;
        QChar postResidue;

        QHash<ResidueIndex, ModificationMass> modifications;

        static Err modificationFromString(
                const QString &modificationString,
                QHash<ResidueIndex, ModificationMass> *result
                ) {

            ERR_INIT

            const QStringList individualModStrings
                    = modificationString.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

            for (const QString &mod : individualModStrings) {

                const QStringList modSplit
                        = mod.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP, QString::SkipEmptyParts);

                e = ErrorUtils::isEqual(modSplit.size(), 2); ree;

                ResidueIndex residueIndex;
                e = ErrorUtils::toInt(modSplit.front(), &residueIndex); ree;

                ModificationMass modificationMass;
                e = ErrorUtils::toDouble(modSplit.back(), &modificationMass); ree;

                result->insert(residueIndex, modificationMass);

            }

            ERR_RETURN
        }

        static QString peptideStringWithMods(
                const QString &sequence,
                QHash<ResidueIndex, ModificationMass> mods
                ) {

            QString peptideBuilder;
            for (int i = 0; i < sequence.size(); i++) {

                const QChar &aa = sequence.at(i);

                peptideBuilder += aa;

                if (mods.contains(i)) {
                    const QString modAddition = QStringLiteral("[%1]").arg(mods.value(i));
                    peptideBuilder += modAddition;
                }
            }

            return peptideBuilder;
        }

        static Err decodePeptideStringWithMods(
                const QString &peptideStringWithMods,
                QString *peptideSequence,
                QString *modString
                ) {

            ERR_INIT

            bool bracketStart = false;
            QString modMassBuilder;

            QStringList mods;

            for (const QChar &c : peptideStringWithMods) {

                if (c == '[') {
                    bracketStart = true;
                    modMassBuilder.clear();
                }

                else if (c == ']') {

                    double modMassString;
                    e = ErrorUtils::toDouble(modMassBuilder, &modMassString); ree;

                    const int modIndex = peptideSequence->size() - 1;

                    const QString mod = S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT
                            .arg(modIndex)
                            .arg(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP)
                            .arg(modMassString);

                    mods.push_back(mod);
                    bracketStart = false;
                }

                else if (bracketStart) {
                    modMassBuilder += c;
                }

                else {
                    peptideSequence->push_back(c);
                }

            }

            *modString = mods.join(S_GLOBAL_SETTINGS.SEPARATOR);

            ERR_RETURN
        }

        [[nodiscard]] QString modificationString() const {

            QStringList modStringBuilder;
            for (auto it = modifications.begin(); it !=  modifications.end(); it++) {

                const QString modIndexString = QString::number(it.key());
                const QString modMassString = QString::number(it.value());
                const QString mod = S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT.arg(
                        modIndexString,
                        S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                        modMassString
                );
                modStringBuilder.push_back(mod);
            }

            return modStringBuilder.join(S_GLOBAL_SETTINGS.SEPARATOR);
        }

        [[nodiscard]] QString peptideStringWithMods() const {
            return peptideStringWithMods(sequence, modifications);
        }

    friend QDataStream &operator <<(QDataStream &stream, const Peptide &pep) {
        stream << pep.id;
        stream << pep.sequence;
        stream << pep.mass;
        stream << pep.isDecoy;
        stream << pep.previousResidue;
        stream << pep.postResidue;
        stream << pep.modifications;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, Peptide &pep) {
        stream >> pep.id;
        stream >> pep.sequence;
        stream >> pep.mass;
        stream >> pep.isDecoy;
        stream >> pep.previousResidue;
        stream >> pep.postResidue;
        stream >> pep.modifications;

        return stream;
    }

};


class ALGORITHMSLIB_EXPORTS PeptidesLibraryTron {

public:

    friend class PeptidesLibraryTronTests;

    explicit PeptidesLibraryTron(const PythiaParameters &pythiaParameters);

    PeptidesLibraryTron() = default;
    ~PeptidesLibraryTron() = default;

    Err exec(
            const QString &fastaFileUri,
            int seed
            );

    [[nodiscard]] QVector<Peptide> peptides() const;

    static Err writePeptidesLib(
            const QVector<Peptide> &peptides,
            const QString &peptidesLibFilePath
            );

    Err readPeptidesLib(const QString &peptidesLibFilePath);

    Err getPeptideById(
            PeptideId peptideId,
            Peptide *peptide
            );

private:

    Err processFasta(const QString &fastaFileUri);
    Err removeDuplicatesFromPeptide();
    Err addDecoys(int seed);
    Err buildPeptides();
    Err addVariableModificationsToPeptides();
    Err addTerminalModificationsToPeptides();
    Err addPeptideIdToPeptides();
    Err addMassToPeptides();
    Err buildPeptidesLookupByPeptideId();

private:

    PythiaParameters m_pythiaParameters;
    QVector<PeptideSequence> m_peptideSequences;
    QVector<Peptide> m_peptides;
    QMap<PeptideId, Peptide> m_peptidesLookupByPeptideId;

};


#endif //PYTHIADIACPP_PEPTIDESLIBRARYTRON_H
