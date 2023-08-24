//
// Created by anichols on 8/23/23.
//

#include "LibraryReader.h"

#include "LibraryCommon.h"

LibraryReader::LibraryReader()
: m_iRT_min(0.0)
, m_iRT_max(0.0)
{}

bool LibraryReader::load(const char *file_name) {

    m_name = std::string(file_name);

    if (LibraryCommon::get_extension(m_name) == std::string(".speclib")) {

        std::ifstream speclib(file_name, std::ifstream::binary);

        if (speclib.fail()) {
            qDebug() << "cannot read the file\n";
            return false;
        }
        read(speclib);
        speclib.close();
    }

    return true;
}

void LibraryReader::read(std::ifstream &in) {

    int gd = 0;
    int gc = 0;
    int ip = 0;
    int version = 0;

    in.read((char*)&version, sizeof(int));

    if (version >= 0) {
        gd = version;
    }
    else {
        in.read((char*)&gd, sizeof(int));
    }

    in.read((char*)&gc, sizeof(int));
    in.read((char*)&ip, sizeof(int));

    ReadWriteCommonFunctons::read_string(in, m_name);
    ReadWriteCommonFunctons::read_string(in, m_fasta_names);
    ReadWriteCommonFunctons::read_array(in, m_proteins);
    ReadWriteCommonFunctons::read_array(in, m_protein_ids);
    ReadWriteCommonFunctons::read_strings(in, m_precursors);
    ReadWriteCommonFunctons::read_strings(in, m_names);
    ReadWriteCommonFunctons::read_strings(in, m_genes);
    in.read((char*)&m_iRT_min, sizeof(double));
    in.read((char*)&m_iRT_max, sizeof(double));
    ReadWriteCommonFunctons::read_array(in, m_entries);

    if (version <= -1 && in.peek() != std::char_traits<char>::eof()) {
        ReadWriteCommonFunctons::read_vector(in, m_elution_groups);
    }
}

const std::vector<Entry> &LibraryReader::getEntries() const {
    return m_entries;
}


