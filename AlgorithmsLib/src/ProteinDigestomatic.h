//
// Created by Drucifer on 11/20/2021.
//

#ifndef PYTHIACPP_PROTEINDIGESTOMATIC_H
#define PYTHIACPP_PROTEINDIGESTOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"
#include "PythiaParameterReader.h"

#include <QStringList>


using namespace Error;

namespace PeptideSequenceNamespace {

    const QString SEQUENCE = QStringLiteral("sequence");
    const QString PREV_RES = QStringLiteral("previousResidue");
    const QString FIRST_RES = QStringLiteral("firstResidue");
    const QString LAST_RES = QStringLiteral("lastResidue");
    const QString POST_RES = QStringLiteral("postResidue");
    const QString MASS = QStringLiteral("mass");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString START_IND = QStringLiteral("startIndex");
    const QString END_IND = QStringLiteral("endIndex");
    const QString SIZE = QStringLiteral("size");

    const QStringList keysToCheck = {
            SEQUENCE,
            PREV_RES,
            FIRST_RES,
            LAST_RES,
            POST_RES,
            MASS,
            IS_DECOY,
            START_IND,
            END_IND,
            SIZE
    };
}


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
            {PeptideSequenceNamespace::SEQUENCE, QVariant(sequence)},
            {PeptideSequenceNamespace::PREV_RES, QVariant(previousResidue)},
            {PeptideSequenceNamespace::FIRST_RES, QVariant(firstResidue)},
            {PeptideSequenceNamespace::LAST_RES, QVariant(lastResidue)},
            {PeptideSequenceNamespace::POST_RES, QVariant(postResidue)},
            {PeptideSequenceNamespace::MASS, QVariant(mass)},
            {PeptideSequenceNamespace::IS_DECOY, QVariant(isDecoy)},
            {PeptideSequenceNamespace::START_IND, QVariant(startIndex)},
            {PeptideSequenceNamespace::END_IND, QVariant(endIndex)},
            {PeptideSequenceNamespace::SIZE, QVariant(size())}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const QStringList &mapKeys = dataMap.keys();
        auto keyCheckLogic = [mapKeys](const QString &s){return mapKeys.contains(s);};
        const bool allKeysPresent = std::all_of(
                PeptideSequenceNamespace::keysToCheck.begin(),
                PeptideSequenceNamespace::keysToCheck.end(),
                keyCheckLogic
                );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        sequence = dataMap.value(PeptideSequenceNamespace::SEQUENCE).toString();
        previousResidue = dataMap.value(PeptideSequenceNamespace::PREV_RES).toChar();
        firstResidue = dataMap.value(PeptideSequenceNamespace::FIRST_RES).toChar();
        lastResidue = dataMap.value(PeptideSequenceNamespace::LAST_RES).toChar();
        postResidue = dataMap.value(PeptideSequenceNamespace::POST_RES).toChar();
        mass = dataMap.value(PeptideSequenceNamespace::MASS).toDouble();
        isDecoy = dataMap.value(PeptideSequenceNamespace::IS_DECOY).toBool();
        startIndex = dataMap.value(PeptideSequenceNamespace::START_IND).toInt();
        endIndex = dataMap.value(PeptideSequenceNamespace::END_IND).toInt();

        ERR_RETURN
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
