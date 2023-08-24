//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_ISOFORM_H
#define PYTHIADIACPP_ISOFORM_H

#include "KarnnLib_Exports.h"

#include "ReadWriteCommonFunctons.h"

#include <fstream>
#include <set>
#include <string>


class KARNNLIB_EXPORTS Isoform {

public:

    Isoform();
    ~Isoform() = default;

    explicit Isoform(const std::string &id);

    Isoform(
        const std::string &id,
        const std::string &name,
        const std::string &gene,
        const std::string &description,
        bool swissprot
        );

    [[nodiscard]] std::string getId() const;

    [[nodiscard]] int getNameIndex() const;

    [[nodiscard]] int getGeneIndex() const;

    friend inline bool operator < (const Isoform &l, const Isoform &r) { return l.getId() < r.getId(); }

    void read(std::ifstream &in);

private:

    std::string m_id;

private:
    std::string m_name;
    std::string m_gene;
    std::string m_description;
    std::set<int> m_precursors;
    int m_nameIndex;
    int m_geneIndex;
    bool m_swissprot;

};


#endif //PYTHIADIACPP_ISOFORM_H
