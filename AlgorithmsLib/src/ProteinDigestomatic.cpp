//
// Created by Drucifer on 11/20/2021.
//

#include "ProteinDigestomatic.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "MathUtils.h"


ProteinDigestomatic::ProteinDigestomatic(const PythiaParameters &params)
: m_pythiaParams(params)
{}


namespace {

    const QChar terminalResidue = '*';

    void setPeptideSequenceAux(
            int startIndex,
            QChar aaPrev,
            QChar aaPost,
            PeptideSequence *peptideSequence
            ) {

        peptideSequence->startIndex = startIndex;
        peptideSequence->endIndex = startIndex + peptideSequence->size() -1;

        peptideSequence->previousResidue = aaPrev;
        peptideSequence->firstResidue = peptideSequence->sequence.front();
        peptideSequence->lastResidue = peptideSequence->sequence.back();
        peptideSequence->postResidue = aaPost;
    }

    void filterSequenceLength(
            int minSequenceLength,
            int maxSequenceLength,
            QVector<PeptideSequence> *peptideSequences
            ) {

        const auto terminatorLogic = [minSequenceLength, maxSequenceLength](const PeptideSequence &seq){
            return seq.size() < minSequenceLength || seq.size() > maxSequenceLength;
        };

        const auto terminator = std::remove_if(
                peptideSequences->begin(),
                peptideSequences->end(),
                terminatorLogic
                );

        peptideSequences->erase(terminator, peptideSequences->end());
    }

    void filterInvalidSequences(QVector<PeptideSequence> *peptideSequences) {

        const QVector<QChar> invalidAminoAcidLetters = {'B', 'J', 'O', 'U', 'Z'};

        const auto charInStringLogic = [invalidAminoAcidLetters](const QChar &c){
            return invalidAminoAcidLetters.contains(c);
        };

        const auto terminatorLogic = [charInStringLogic](const PeptideSequence &ps){
            return std::any_of(
                    ps.sequence.begin(),
                    ps.sequence.end(),
                    charInStringLogic
                    );
        };

        const auto terminator = std::remove_if(
                peptideSequences->begin(),
                peptideSequences->end(),
                terminatorLogic
                );

        peptideSequences->erase(terminator, peptideSequences->end());
    }


}//END NAMESPACE
Err ProteinDigestomatic::digestProtein(
        const ProteinSequence &_proteinSequence,
        QVector<PeptideSequence> *peptideSequences
        ) {

    ERR_INIT;

    e = ErrorUtils::isNotEmpty(_proteinSequence); ree;

    ProteinSequence proteinSequence = _proteinSequence;

    if (m_pythiaParams.replaceLeucinesWithX) {
        proteinSequence = proteinSequence.replace('L', 'X');
        proteinSequence = proteinSequence.replace('I', 'X');
    }

    PeptideSequence peptideSequence;

    int startIndex = 0;
    for(int i = 0; i < proteinSequence.size(); ++i){

        const QChar &aa = proteinSequence.at(i);

        if(m_pythiaParams.nTermCleavePoints.contains(aa)){

            const QChar &aaPrev =  startIndex == 0 ? terminalResidue : proteinSequence.at(startIndex-1);
            const QChar &aaPost = aa;

            if(!peptideSequence.sequence.isEmpty()){

                setPeptideSequenceAux(
                        startIndex,
                        aaPrev,
                        aaPost,
                        &peptideSequence
                        );

                peptideSequences->push_back(peptideSequence);
            }

            peptideSequence = PeptideSequence();
            startIndex = i;

            peptideSequence.sequence += aa.toLatin1();
        }

        else if (m_pythiaParams.cTermCleavePoints.contains(aa) && i < proteinSequence.size() - 1){

            const QChar &aaPrev =  startIndex == 0 ? terminalResidue : proteinSequence.at(startIndex-1);
            const QChar &aaPost = proteinSequence.at(i+1);

            peptideSequence.sequence += aa.toLatin1();

            setPeptideSequenceAux(
                    startIndex,
                    aaPrev,
                    aaPost,
                    &peptideSequence
                    );

            peptideSequences->push_back(peptideSequence);

            peptideSequence = PeptideSequence();
            startIndex = i + 1;
        }

        else{
            peptideSequence.sequence += aa.toLatin1();
        }

    }

    //Handle last peptide edge case.
    if(!peptideSequence.sequence.isEmpty()){

        const int startIndexCheck = startIndex-1;
        const QChar &aaPrev = startIndexCheck < 0 ? terminalResidue : proteinSequence.at(startIndexCheck);
        const QChar &aaPost = terminalResidue;
        setPeptideSequenceAux(startIndex,
                              aaPrev,
                              aaPost,
                              &peptideSequence);

        peptideSequences->push_back(peptideSequence);
    }

    e = createPartialDigestPeptides(peptideSequences); ree;
    e = createRaggedSegments(peptideSequences); ree;

    filterSequenceLength(m_pythiaParams.peptideLengthMin,
                         m_pythiaParams.peptideLengthMax,
                         peptideSequences);

    filterInvalidSequences(peptideSequences);
    calculateMasses(peptideSequences);

    ERR_RETURN;
}

Err ProteinDigestomatic::createPartialDigestPeptides(QVector<PeptideSequence> *peptideSequences) const
{
    ERR_INIT;

    e = ErrorUtils::isNotEmpty(*peptideSequences); ree;

    if (m_pythiaParams.allowedMissedCleavages < 1){
        ERR_RETURN;
    }

    QVector<PeptideSequence> peptideSequencePartials;

    for(int i = 1; i <= m_pythiaParams.allowedMissedCleavages; ++i){

        const int firstArraySliceSize = peptideSequences->size() - i;
        if(firstArraySliceSize < 1){
            continue;
        }

        QVector<PeptideSequence> array1 = peptideSequences->mid(0, firstArraySliceSize);

        for(int j = 1; j < i + 1; j++){

            const QVector<PeptideSequence> arrayToAdd = peptideSequences->mid(j, firstArraySliceSize);

            if(array1.size() != arrayToAdd.size()){
                qDebug() << "SISD" << peptideSequences->size();
                qDebug() << 0 << firstArraySliceSize;
                qDebug() << j << firstArraySliceSize;
            }

            e = ErrorUtils::isEqual(array1.size(), arrayToAdd.size()); ree;

            for(int k = 0; k < array1.size(); k++){
                PeptideSequence &a1 = array1[k];
                a1.sequence += arrayToAdd.at(k).sequence;
                a1.endIndex = arrayToAdd.at(k).endIndex;
                a1.lastResidue = arrayToAdd.at(k).lastResidue;
                a1.postResidue = arrayToAdd.at(k).postResidue;
            }
        }

        peptideSequencePartials.append(array1);
    }

    peptideSequences->append(peptideSequencePartials);
    ERR_RETURN;
}

Err ProteinDigestomatic::createRaggedSegments(QVector<PeptideSequence> *peptideSequences) const {

    ERR_INIT;

    e = ErrorUtils::isNotEmpty(*peptideSequences); ree;

    QVector<PeptideSequence> tempPeptideSequence;
    QSet<QString> raggedAlreadyExists;

    for(const PeptideSequence &pepSeq : *peptideSequences){

        const int seqSize = pepSeq.sequence.size() - 2;
        const QStringList seqSplit = pepSeq.sequence.split("", QString::SkipEmptyParts);

        QStringList cRaggedBuilder = seqSplit;
        QStringList nRaggedBuilder = seqSplit;

        //-1 is used here so we don't add single residues as segments.
        for(int i = 0; i < seqSize; ++i){

            if(m_pythiaParams.raggedness == Raggedness::BothRagged
               || m_pythiaParams.raggedness == Raggedness::CTermRagged){

                const QChar postResidue = cRaggedBuilder.back()[0];
                cRaggedBuilder.pop_back();
                const QString poppedSeq = cRaggedBuilder.join("");

                PeptideSequence newPepSeq;
                newPepSeq.startIndex = pepSeq.startIndex;
                newPepSeq.endIndex = pepSeq.endIndex - i - 1;
                newPepSeq.sequence = poppedSeq;
                newPepSeq.firstResidue = newPepSeq.sequence.front();
                newPepSeq.lastResidue = newPepSeq.sequence.back();
                newPepSeq.previousResidue = pepSeq.previousResidue;
                newPepSeq.postResidue = postResidue;
                newPepSeq.isDecoy = pepSeq.isDecoy;

                if(!newPepSeq.sequence.isEmpty()
                   && !raggedAlreadyExists.contains(newPepSeq.sequence)){
                    tempPeptideSequence.push_back(newPepSeq);
                    raggedAlreadyExists.insert(newPepSeq.sequence);
                }
            }

            if(m_pythiaParams.raggedness == Raggedness::BothRagged
               || m_pythiaParams.raggedness == Raggedness::NTermRagged){

                const QChar preResidue = nRaggedBuilder.front()[0];
                nRaggedBuilder.pop_front();
                const QString poppedSeq = nRaggedBuilder.join("");

                PeptideSequence newPepSeq;
                newPepSeq.startIndex = pepSeq.startIndex + i + 1;
                newPepSeq.endIndex = pepSeq.endIndex;
                newPepSeq.sequence = poppedSeq;
                newPepSeq.firstResidue = newPepSeq.sequence.front();
                newPepSeq.lastResidue = newPepSeq.sequence.back();
                newPepSeq.previousResidue =preResidue;
                newPepSeq.postResidue = pepSeq.postResidue;
                newPepSeq.isDecoy = pepSeq.isDecoy;

                if(!newPepSeq.sequence.isEmpty()
                   && !raggedAlreadyExists.contains(newPepSeq.sequence)){
                    tempPeptideSequence.push_back(newPepSeq);
                    raggedAlreadyExists.insert(newPepSeq.sequence);
                }
            }
        }
    }

    peptideSequences->append(tempPeptideSequence);
    ERR_RETURN;
}

void ProteinDigestomatic::calculateMasses(QVector<PeptideSequence> *peptideSequences) {

    for (int i = 0; i < peptideSequences->size(); i++) {

        PeptideSequence &ps = (*peptideSequences)[i];
        ps.mass = MathUtils::pRound(BiophysicalCalcs::calculatePeptideMass(
                ps.sequence,
                m_pythiaParams.aminoAcids
                ), 4);
    }

}
