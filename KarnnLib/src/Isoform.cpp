//
// Created by anichols on 8/23/23.
//

#include "Isoform.h"

Isoform::Isoform()
: m_nameIndex(0)
, m_geneIndex(0)
, m_swissprot(true)
{}

Isoform::Isoform(const std::string &id)
: m_id(id)
, m_nameIndex(0)
, m_geneIndex(0)
, m_swissprot(true)
{}

Isoform::Isoform(
        const std::string &id,
        const std::string &name,
        const std::string &gene,
        const std::string &description,
        bool swissprot
        )
: m_id(id)
, m_name(name)
, m_gene(gene)
, m_description(description)
, m_swissprot(swissprot)
, m_nameIndex(0)
, m_geneIndex(0)
{}

void Isoform::read(std::ifstream &in) {

    int sp = 0, size = 0;
    in.read((char*)&sp, sizeof(int));
    in.read((char*)&size, sizeof(int));
    m_swissprot = sp;

    ReadWriteCommonFunctons::read_string(in, m_id);
    ReadWriteCommonFunctons::read_string(in, m_name);
    ReadWriteCommonFunctons::read_string(in, m_gene);
    in.read((char*)&m_nameIndex, sizeof(int));
    in.read((char*)&m_geneIndex, sizeof(int));
    m_precursors.clear();
    for (int i = 0; i < size; i++) {

        int pr = -1;

        in.read((char*)&pr, sizeof(int));

        if (pr >= 0) {
            m_precursors.insert(pr);
        }
    }
}

std::string Isoform::getId() const {
    return m_id;
}

int Isoform::getNameIndex() const {
    return m_nameIndex;
}

int Isoform::getGeneIndex() const {
    return m_geneIndex;
}

