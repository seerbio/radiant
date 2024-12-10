//
// Created by anichols on 6/1/24.
//

#include "SequenceSubstringSearchomatic.h"

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>


typedef seqan::String<seqan::AminoAcid> TSequence;
typedef seqan::StringSet<TSequence> TStringSet;
typedef seqan::Index<TStringSet, seqan::IndexEsa<> > TIndex;


Err SequenceSubstringSearchomatic::init(const QVector<FastaEntry> &fastaEntries) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fastaEntries); ree;
    m_fastaEntries = fastaEntries;

    ERR_RETURN
}

Err SequenceSubstringSearchomatic::findProteinGroups(
    const QStringList &searchSequences,
    QVector<QString> *proteinGroupsAll
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(searchSequences); ree;

    TStringSet sequences;
    for (const FastaEntry &fe : m_fastaEntries) {
        QString seq = fe.fastaSequence;
        seq = seq.replace("I", "L");
        appendValue(sequences, seqan::CharString(seq.toStdString()));
    }

    TIndex index(sequences);

    QVector<QVector<int>> proteinDescriptionIndexesAll;
    for (const QString &searchSequence : searchSequences) {

        QString searchSequenceIsoLeucinesReplace = searchSequence;
        searchSequenceIsoLeucinesReplace = searchSequenceIsoLeucinesReplace.replace("I", "L");

        const TSequence pattern = searchSequenceIsoLeucinesReplace.toStdString();

        seqan::Finder<TIndex> finder(index);

        QVector<int> proteinDescriptionIndexes;
        while (find(finder, pattern)) {

            auto pos = seqan::position(finder);
            unsigned seqNo = seqan::getSeqNo(pos);
            unsigned seqOffset = seqan::getSeqOffset(pos);

            if (seqan::infix(sequences[seqNo], seqOffset, seqOffset + length(pattern)) == pattern) {
                proteinDescriptionIndexes.push_back(static_cast<int>(seqNo));
            }
        }

        proteinDescriptionIndexesAll.push_back(proteinDescriptionIndexes);
    }

    e = ErrorUtils::isEqual(searchSequences.size(), proteinDescriptionIndexesAll.size()); ree;

    for (int i = 0; i < searchSequences.size(); i++) {

        const QVector<int> &proteinDescriptionIndexes = proteinDescriptionIndexesAll.at(i);

        if (proteinDescriptionIndexes.isEmpty()) {
            qWarning() << searchSequences.at(i) << "not found in fasta for reannotation";
        }

        QStringList proteinGroups;
        for (int ind : proteinDescriptionIndexes) {
            proteinGroups.push_back(m_fastaEntries.at(ind).fastaDescription);
        }

        QString proteinGroupsString = proteinGroups.join(';');
        proteinGroupsString = proteinGroupsString.replace(",", " ");
        proteinGroupsAll->push_back(proteinGroupsString);
    }

    ERR_RETURN
}
