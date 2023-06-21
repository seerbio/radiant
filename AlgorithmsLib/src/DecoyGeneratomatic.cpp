//
// Created by anichols on 6/21/23.
//

#include "DecoyGeneratomatic.h"

#include "ErrorUtils.h"

#include <random>


namespace {

    QString shufflePeptide(const QString& peptideSeq) {

        std::mt19937 rng(666); //TODO make this settable.

        const int minSizePeptideShuffle = 5;
        const int peptideLen = peptideSeq.size();
        if (peptideLen < minSizePeptideShuffle) {
            return peptideSeq;
        }

        QStringList strList = peptideSeq.split("", QString::SkipEmptyParts);

        std::reverse(strList.rbegin() + 1, strList.rend() - 1);

        std::shuffle(
                strList.begin() + 1,
                strList.begin() + peptideLen - 1,
                rng
        );

        return strList.join("");
    }

}//namespace
QString DecoyGeneratomatic::midShuffleLogicDecoysLogic(const QString &peptideSeq) {

    QString newPepSeq = peptideSeq;

    int retries = 0;
    const int maxRetryCount = 10;
    while (newPepSeq == peptideSeq) {

        newPepSeq = shufflePeptide(peptideSeq);

        if (retries++ > maxRetryCount) {
            newPepSeq = S_GLOBAL_SETTINGS.FAILED_SHUFFLE;
            return newPepSeq;
        }
    }

    return newPepSeq;
}

QString DecoyGeneratomatic::diannDecoyLogic(const QString &peptideSeq) {

    const int minPeptideSizeForDecoyification = 4;
    if (peptideSeq.size() < minPeptideSizeForDecoyification) {
        return peptideSeq;
    }

    const QMap<QChar, QChar> aminoAcidSubstitutionMap = {
            {'G', 'L'},
            {'A', 'L'},
            {'V', 'L'},
            {'L', 'V'},
            {'I', 'V'},
            {'X', 'V'},
            {'J', 'V'},
            {'F', 'L'},
            {'M', 'L'},
            {'P', 'L'},
            {'W', 'L'},
            {'S', 'T'},
            {'C', 'S'},
            {'T', 'S'},
            {'Y', 'S'},
            {'H', 'S'},
            {'K', 'L'},
            {'R', 'L'},
            {'Q', 'N'},
            {'E', 'D'},
            {'N', 'E'},
            {'D', 'Q'}
    };

    QString newPeptideSequence = peptideSeq;
    const int indexCTermSub = peptideSeq.size() - 2;

    newPeptideSequence[1] = aminoAcidSubstitutionMap[peptideSeq[1]];
    newPeptideSequence[indexCTermSub] = aminoAcidSubstitutionMap[peptideSeq[indexCTermSub]];

    return newPeptideSequence;
}
