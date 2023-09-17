//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_LIBRARYREADERR_H
#define PYTHIADIACPP_LIBRARYREADER_H

#include "KarnnLib_Exports.h"

#include "Entry.h"
#include "Isoform.h"
#include "ProteinGroup.h"
#include "ReadWriteCommonFunctons.h"

#include <QDebug>

#include <fstream>
#include <set>
#include <string>
#include <vector>


class KARNNLIB_EXPORTS LibraryReader {

    using PG = ProteinGroup;

public:

    LibraryReader();
    ~LibraryReader() = default;

    bool load(const std::string &fileName);

    [[nodiscard]] const std::vector<Entry> &getEntries() const;


private:

    void read(std::ifstream &in);

private:
    std::string m_name;
    std::string m_fasta_names;
    std::vector<Isoform> m_proteins;
    std::vector<PG> m_protein_ids;
    std::vector<PG> m_protein_groups;
    std::vector<PG> m_gene_groups;
    std::vector<int> m_gg_index;
    std::vector<std::string> m_precursors;
    std::vector<std::string> m_names;
    std::vector<std::string> m_genes;

    double m_iRT_min;
    double m_iRT_max;

    std::map<std::string, int> m_eg;
    std::vector<int> m_elution_groups;
    std::vector<int> m_co_elution;
    std::vector<std::pair<int, int> > m_co_elution_index;

    std::vector<Entry> m_entries;

};


#endif //PYTHIADIACPP_LIBRARYREADERR_H
