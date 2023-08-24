//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_PROTEINGROUP_H
#define PYTHIADIACPP_PROTEINGROUP_H


#include "KarnnLib_Exports.h"

#include "ReadWriteCommonFunctons.h"

#include <set>
#include <string>
#include <vector>

class Isoform;

class KARNNLIB_EXPORTS ProteinGroup {

public:

    using PG = ProteinGroup;

    std::string ids;
    mutable std::string names;
    mutable std::string genes;
    mutable std::vector<int> precursors; // precursor indices in the library
    mutable std::set<int> proteins;
    std::vector<int> name_indices;
    std::vector<int> gene_indices;

    ProteinGroup() = default;

    ~ProteinGroup() = default;

    explicit ProteinGroup(std::string ids);

    friend inline bool operator < (const PG &left, const PG &right) { return left.ids < right.ids; }

    void annotate(
            std::vector<Isoform> &proteins,
            std::vector<std::string> &names,
            std::vector<std::string> &genes
            );

    template <class F>
    void write(F &out) {
        int size_p = proteins.size();
        out.write((char*)&size_p, sizeof(int));
        ReadWriteCommonFunctons::write_string(out, ids);
        ReadWriteCommonFunctons::write_string(out, names);
        ReadWriteCommonFunctons::write_string(out, genes);
        ReadWriteCommonFunctons::write_vector(out, name_indices);
        ReadWriteCommonFunctons::write_vector(out, gene_indices);
        ReadWriteCommonFunctons::write_vector(out, precursors);

        for (auto it = proteins.begin(); it != proteins.end(); it++) {
            out.write((char*)&(*it), sizeof(int));
        }
    }

    template <class F>
    void read(F &in) {
        int sp = 0, size_p = 0;
        in.read((char*)&size_p, sizeof(int));

        ReadWriteCommonFunctons::read_string(in, ids);
        ReadWriteCommonFunctons::read_string(in, names);
        ReadWriteCommonFunctons::read_string(in, genes);
        ReadWriteCommonFunctons::read_vector(in, name_indices);
        ReadWriteCommonFunctons::read_vector(in, gene_indices);
        ReadWriteCommonFunctons::read_vector(in, precursors);
        for (int i = 0; i < size_p; i++) {
            int p = -1;
            in.read((char*)&p, sizeof(int));
            if (p >= 0) proteins.insert(p);
        }
    }
};

#endif //PYTHIADIACPP_PROTEINGROUP_H
