//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_ISOFORM_H
#define PYTHIADIACPP_ISOFORM_H

#include "KarnnLib_Exports.h"

#include "ReadWriteCommonFunctons.h"

#include <set>
#include <string>


class KARNNLIB_EXPORTS Isoform {

public:

    std::string id;
    mutable std::string name;
    mutable std::string gene;
    mutable std::string description;
    mutable std::set<int> precursors; // precursor indices in the library
    int name_index = 0;
    int gene_index = 0;
    bool swissprot = true;

    Isoform() = default;
    Isoform(const std::string &id);
    Isoform(
            const std::string &id,
            const std::string &name,
            const std::string &gene,
            const std::string &description,
            bool swissprot
            );

    friend inline bool operator < (const Isoform &left, const Isoform &right) { return left.id < right.id; }

    template <class F>
    void write(F &out) {
        int sp = swissprot;
        int size = static_cast<int>(precursors.size());
        out.write((char*)&sp, sizeof(int));
        out.write((char*)&size, sizeof(int));
        ReadWriteCommonFunctons::write_string(out, id);
        ReadWriteCommonFunctons::write_string(out, name);
        ReadWriteCommonFunctons::write_string(out, gene);
        out.write((char*)&name_index, sizeof(int));
        out.write((char*)&gene_index, sizeof(int));

        for (auto it = precursors.begin(); it != precursors.end(); it++) {
            out.write((char*)&(*it), sizeof(int));
        }
    }

    template <class F>
    void read(F &in) {
        int sp = 0, size = 0;
        in.read((char*)&sp, sizeof(int));
        in.read((char*)&size, sizeof(int));
        swissprot = sp;

        ReadWriteCommonFunctons::read_string(in, id);
        ReadWriteCommonFunctons::read_string(in, name);
        ReadWriteCommonFunctons::read_string(in, gene);
        in.read((char*)&name_index, sizeof(int));
        in.read((char*)&gene_index, sizeof(int));
        precursors.clear();
        for (int i = 0; i < size; i++) {
            int pr = -1;
            in.read((char*)&pr, sizeof(int));
            if (pr >= 0) precursors.insert(pr);
        }
    }
};


#endif //PYTHIADIACPP_ISOFORM_H
