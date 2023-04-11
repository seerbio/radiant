//
// Created by anichols on 4/10/23.
//

#ifndef PYTHIADIACPP_PEPTIDEMASSRTREE_H
#define PYTHIADIACPP_PEPTIDEMASSRTREE_H


#include "AlgorithmsLib_Exports.h"
#include "AminoAcids.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;


class ALGORITHMSLIB_EXPORTS PeptideMassRTree {

public:

    PeptideMassRTree();
    ~PeptideMassRTree();

    Err init(
            const QVector<PeptideStringWithMods> &peptideStringsWithMods,
            const AminoAcids &aminoAcids
            );

private:

    Q_DISABLE_COPY(PeptideMassRTree) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_PEPTIDEMASSRTREE_H
