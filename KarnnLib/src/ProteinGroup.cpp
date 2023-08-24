//
// Created by anichols on 8/23/23.
//

#include "ProteinGroup.h"

#include "Isoform.h"

#include "cpu_info.h"


ProteinGroup::ProteinGroup(std::string ids) : ids(ids) {}

void ProteinGroup::annotate(
        std::vector<Isoform> &_proteins,
        std::vector<std::string> &_names,
        std::vector<std::string> &_genes
        ) {

    std::set<int> name;
    std::set<int> gene;

    std::string word;

    if (!_names.empty()) {

        for (auto &p : proteins) {
            if (!_names[_proteins[p].getNameIndex()].empty()) {
                name.insert(_proteins[p].getNameIndex());
            }
        }

        if (!names.empty() && proteins.empty()) {

            std::stringstream list(names);
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
        names.clear(); name_indices.clear(); name_indices.insert(name_indices.begin(), name.begin(), name.end());
        if (!name_indices.empty()) {
            for (int i = 0; i < name_indices.size(); i++) {
                names += (names.size() ? std::string(";") : std::string("")) + _names[name_indices[i]];
            }
        }
    }
    if (!_genes.empty()) {
        for (auto &p : proteins) {
            if (!_genes[_proteins[p].getGeneIndex()].empty()) {
                gene.insert(_proteins[p].getGeneIndex());
            }
        }
        if (genes.size() && !proteins.size()) {
            std::stringstream list(genes);
            while (std::getline(list, word, ';')) {
                word = fast_trim(word);
                auto pos = std::lower_bound(_genes.begin(), _genes.end(), word);
                if (pos != _genes.end()) if (*pos == word) gene.insert(std::distance(_genes.begin(), pos));
            }
        }
        genes.clear(); gene_indices.clear(); gene_indices.insert(gene_indices.begin(), gene.begin(), gene.end());
        if (gene_indices.size()) for (int i = 0; i < gene_indices.size(); i++) genes += (genes.size() ? std::string(";") : std::string("")) + _genes[gene_indices[i]];
    }
}
