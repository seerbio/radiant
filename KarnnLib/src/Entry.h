//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_ENTRY_H
#define PYTHIADIACPP_ENTRY_H


#include "Peptide.h"
#include "ProteinGroup.h"

#include <set>


class Entry {

    using PG = ProteinGroup;

public:

    Peptide target;
    Peptide decoy;
    std::string name; // precursor id
    std::set<PG>::iterator prot;
    int pid_index = 0;

    friend bool operator < (const Entry &left, const Entry &right) { return left.name < right.name; }

    template <class F>
    void read(F &in) {
        target.read(in);
        int dc = 0; in.read((char*)&dc, sizeof(int));
        if (dc) decoy.read(in);

        int ff = 0, prt = 0;
        in.read((char*)&ff, sizeof(int));
        in.read((char*)&prt, sizeof(int));
        in.read((char*)&pid_index, sizeof(int));
        ReadWriteCommonFunctons::read_string(in, name);
    }
};


#endif //PYTHIADIACPP_ENTRY_H
