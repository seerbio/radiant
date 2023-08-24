//
// Created by anichols on 8/23/23.
//

#include "ProteinGroup.h"

#include "Isoform.h"

#include "cpu_info.h"


ProteinGroup::ProteinGroup(std::string ids) : m_ids(ids) {}

void ProteinGroup::annotate(
        std::vector<Isoform> &_proteins,
        std::vector<std::string> &_names,
        std::vector<std::string> &_genes
        ) {

    std::set<int> name;
    std::set<int> gene;

    std::string word;

    if (!_names.empty()) {

        for (auto &p : m_proteins) {
            if (!_names[_proteins[p].getNameIndex()].empty()) {
                name.insert(_proteins[p].getNameIndex());
            }
        }

        if (!m_names.empty() && m_proteins.empty()) {

            std::stringstream list(m_names);
            while (std::getline(list, word, ';')) {
                word = fast_trim(word);
                auto pos = std::lower_bound(_names.begin(), _names.end(), word);
                if (pos != _names.end()) {
                    if (*pos == word) {
                        name.insert(std::distance(_names.begin(), pos));
                    }
                }
            }
        }
        m_names.clear(); m_nameIndices.clear(); m_nameIndices.insert(m_nameIndices.begin(), name.begin(), name.end());
        if (!m_nameIndices.empty()) {
            for (int i = 0; i < m_nameIndices.size(); i++) {
                m_names += (m_names.size() ? std::string(";") : std::string("")) + _names[m_nameIndices[i]];
            }
        }
    }

    if (!_genes.empty()) {
        for (auto &p : m_proteins) {
            if (!_genes[_proteins[p].getGeneIndex()].empty()) {
                gene.insert(_proteins[p].getGeneIndex());
            }
        }
        if (!m_genes.empty() && m_proteins.empty()) {
            std::stringstream list(m_genes);
            while (std::getline(list, word, ';')) {
                word = fast_trim(word);
                auto pos = std::lower_bound(_genes.begin(), _genes.end(), word);
                if (pos != _genes.end()) if (*pos == word) gene.insert(std::distance(_genes.begin(), pos));
            }
        }
        m_genes.clear(); m_geneIndices.clear(); m_geneIndices.insert(m_geneIndices.begin(), gene.begin(), gene.end());
        if (!m_geneIndices.empty()) {
            for (int i = 0; i < m_geneIndices.size(); i++) {
                m_genes += (m_genes.size() ? std::string(";") : std::string("")) + _genes[m_geneIndices[i]];
            }
        }
    }
}

void ProteinGroup::read(std::ifstream &in) {

    int sp = 0;
    int size_p = 0;

    in.read((char*)&size_p, sizeof(int));

    ReadWriteCommonFunctons::read_string(in, m_ids);
    ReadWriteCommonFunctons::read_string(in, m_names);
    ReadWriteCommonFunctons::read_string(in, m_genes);
    ReadWriteCommonFunctons::read_vector(in, m_nameIndices);
    ReadWriteCommonFunctons::read_vector(in, m_geneIndices);
    ReadWriteCommonFunctons::read_vector(in, m_precursors);

    for (int i = 0; i < size_p; i++) {
        int p = -1;
        in.read((char*)&p, sizeof(int));
        if (p >= 0) m_proteins.insert(p);
    }
}

std::string ProteinGroup::getIds() const {
    return m_ids;
}

std::string ProteinGroup::getNames() const {
    return m_names;
}

std::string ProteinGroup::getGenes() const {
    return m_genes;
}

std::vector<int> ProteinGroup::getPrecursors() const {
    return m_precursors;
}

std::set<int> ProteinGroup::getProteins() const {
    return m_proteins;
}

std::vector<int> ProteinGroup::getNameIndices() const {
    return m_geneIndices;
}

std::vector<int> ProteinGroup::getGeneIndices() const {
    return m_geneIndices;
}
