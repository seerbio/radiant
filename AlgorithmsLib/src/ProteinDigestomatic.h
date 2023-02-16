//
// Created by Drucifer on 11/20/2021.
//

#ifndef PYTHIACPP_PROTEINDIGESTOMATIC_H
#define PYTHIACPP_PROTEINDIGESTOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "Error.h"

#include <QStringList>


using namespace Error;


struct ALGORITHMSLIB_EXPORTS PeptideSequence {

        QString sequence;
        QChar previousResidue;
        QChar firstResidue;
        QChar lastResidue;
        QChar postResidue;

        bool isDecoy = false;
        int startIndex = -1;
        int endIndex = -1;

        [[nodiscard]] int size() const {
            return sequence.size();
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

private:

    PythiaParameters m_pythiaParams;

};


#endif //PYTHIACPP_PROTEINDIGESTOMATIC_H
