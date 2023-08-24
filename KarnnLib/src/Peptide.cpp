//
// Created by anichols on 8/23/23.
//

#include "Peptide.h"

Peptide::Peptide()
: m_index(0)
, m_charge(0)
, m_length(0)
, m_noCal(0)
, m_mz(0.0)
, m_iRT(0.0)
, m_sRT(0.0)
, m_libQvalue(0.0)
{}

void Peptide::init(
        float _mz,
        float _iRT,
        int _charge,
        int _index
        ) {

    m_mz = _mz;
    m_iRT = _iRT;
    m_sRT = 0.0;
    m_charge = _charge;
    m_index = _index;
}

void Peptide::read(std::ifstream &in) {
    in.read((char*)&m_index, sizeof(int));
    in.read((char*)&m_charge, sizeof(int));
    in.read((char*)&m_length, sizeof(int));
    in.read((char*)&m_mz, sizeof(float));
    in.read((char*)&m_iRT, sizeof(float));
    in.read((char*)&m_sRT, sizeof(float));

    ReadWriteCommonFunctons::read_vector(in, m_fragments);
}

void Peptide::free() {
    std::vector<Product>().swap(m_fragments);
}

int Peptide::getIndex() const {
    return m_index;
}

int Peptide::getCharge() const {
    return m_charge;
}

int Peptide::getLength() const {
    return m_length;
}

int Peptide::getNoCal() const {
    return m_noCal;
}

float Peptide::getMz() const {
    return m_mz;
}

float Peptide::getIRT() const {
    return m_iRT;
}

float Peptide::getSRT() const {
    return m_sRT;
}

float Peptide::getLibQValue() const {
    return m_libQvalue;
}

std::vector<Product> Peptide::getFragments() const {
    return m_fragments;
}
