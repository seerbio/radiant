//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_PROTEINGROUP_H
#define PYTHIADIACPP_PROTEINGROUP_H


#include "KarnnLib_Exports.h"

#include "ReadWriteCommonFunctons.h"

#include <fstream>
#include <set>
#include <string>
#include <vector>

class Isoform;

class KARNNLIB_EXPORTS ProteinGroup {

    using PG = ProteinGroup;

public:

    ProteinGroup() = default;
    ~ProteinGroup() = default;

    explicit ProteinGroup(std::string ids);

    void annotate(
            std::vector<Isoform> &proteins,
            std::vector<std::string> &names,
            std::vector<std::string> &genes
            );

    void read(std::ifstream &in);

    [[nodiscard]] std::string getIds() const;
    [[nodiscard]] std::string getNames() const;
    [[nodiscard]] std::string getGenes() const;
    [[nodiscard]] std::vector<int> getPrecursors() const;
    [[nodiscard]] std::set<int> getProteins() const;
    [[nodiscard]] std::vector<int> getNameIndices() const;
    [[nodiscard]] std::vector<int> getGeneIndices() const;

public:
    friend inline bool operator < (const PG &left, const PG &right) { return left.getIds() < right.getIds();}

private:

    std::string m_ids;
    std::string m_names;
    std::string m_genes;
    std::vector<int> m_precursors;
    std::set<int> m_proteins;
    std::vector<int> m_nameIndices;
    std::vector<int> m_geneIndices;

};

#endif //PYTHIADIACPP_PROTEINGROUP_H
