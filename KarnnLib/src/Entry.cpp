//
// Created by anichols on 8/23/23.
//

#include "Entry.h"

void Entry::read(std::ifstream &in) {

    m_target.read(in);
    int dc = 0; in.read((char*)&dc, sizeof(int));
    if (dc) {
        m_decoy.read(in);
    }

    int ff = 0;
    int prt = 0;
    in.read((char*)&ff, sizeof(int));
    in.read((char*)&prt, sizeof(int));
    in.read((char*)&m_pidIndex, sizeof(int));
    ReadWriteCommonFunctons::read_string(in, m_name);

}

std::string Entry::getName() const {
    return m_name;
}

Peptide Entry::getTarget() const {
    return m_target;
}

Peptide Entry::getDecoy() const {
    return m_decoy;
}
