//
// Created by Drucifer on 11/20/2021.
//

#ifndef PYTHIACPP_PROTEINDIGESTOMATIC_H
#define PYTHIACPP_PROTEINDIGESTOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "ParquetReader.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "Error.h"

#include <QStringList>


using namespace Error;


class ALGORITHMSLIB_EXPORTS PeptideSequence : public ParquetReaderInputBase {

public:
        QString sequence;
        QChar previousResidue;
        QChar firstResidue;
        QChar lastResidue;
        QChar postResidue;

        double mass = -1.0;

        bool isDecoy = false;
        int startIndex = -1;
        int endIndex = -1;

        [[nodiscard]] int size() const {
            return sequence.size();
        }

        QMap<QString, QVariant> map() override {

            return {
                {"sequence", QVariant(sequence)},
                {"previousResidue", QVariant(previousResidue)},
                {"firstResidue", QVariant(firstResidue)},
                {"lastResidue", QVariant(lastResidue)},
                {"postResidue", QVariant(postResidue)},
                {"mass", QVariant(mass)},
//                {"isDecoy", QVariant(isDecoy)},
                {"startIndex", QVariant(startIndex)},
                {"endIndex", QVariant(endIndex)},
                {"size", QVariant(size())}
            };
        }
};


class ALGORITHMSLIB_EXPORTS ProteinDigestomatic {

public:

    explicit ProteinDigestomatic(const PythiaParameters &params);

    Err digestProtein(
            const ProteinSequence &proteinSequence,
            QVector<PeptideSequence> *peptideSequences
            );

private:

    Err createPartialDigestPeptides(QVector<PeptideSequence> *peptideSequences) const;
    Err createRaggedSegments(QVector<PeptideSequence> *peptideSequences) const;
    void calculateMasses(QVector<PeptideSequence> *peptideSequences);

private:

    PythiaParameters m_pythiaParams;

};


#endif //PYTHIACPP_PROTEINDIGESTOMATIC_H
