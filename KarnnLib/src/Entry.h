//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_ENTRY_H
#define PYTHIADIACPP_ENTRY_H

#include "Peptide.h"
#include "ProteinGroup.h"

#include <QDebug>

#include <fstream>
#include <set>


class Entry {

    using PG = ProteinGroup;

public:

    Entry() = default;
    ~Entry() = default;

    [[nodiscard]] std::string getName() const;

    [[nodiscard]] Peptide getTarget() const;

    [[nodiscard]] Peptide getDecoy() const;

    friend bool operator < (const Entry &l, const Entry &r) {return l.getName() < r.getName();}

    void read(std::ifstream &in);

private:

    Peptide m_target;
    Peptide m_decoy;
    std::string m_name;
    std::set<PG>::iterator m_prot;
    int m_pidIndex = 0;

};


#endif //PYTHIADIACPP_ENTRY_H
