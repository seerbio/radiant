//
// Created by anichols on 6/1/24.
//

#ifndef SEQUENCESUBSTRINGSEARCHOMATIC_H
#define SEQUENCESUBSTRINGSEARCHOMATIC_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "FastaReader.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS SequenceSubstringSearchomatic {

public:

    SequenceSubstringSearchomatic() = default;
    ~SequenceSubstringSearchomatic() = default;

    Err init(const QVector<FastaEntry> &fastaEntries);

    Err findProteinGroups(
        const QStringList &searchSequences,
        QVector<QString> *proteinGroupsAll
        );

private:

    QVector<FastaEntry> m_fastaEntries;


};


#endif //SEQUENCESUBSTRINGSEARCHOMATIC_H
