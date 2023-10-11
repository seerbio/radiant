//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_PEPTIDE_H
#define PYTHIADIACPP_PEPTIDE_H

#include "KarnnLib_Exports.h"

#include "Product.h"
#include "ReadWriteCommonFunctons.h"

#include <fstream>
#include <vector>

class KARNNLIB_EXPORTS Peptide {

public:

    Peptide();
    ~Peptide() = default;

    void init(
            float _mz,
            float _iRT,
            int _charge,
            int _index
            );

    void free();

    void read(std::ifstream &in);

    [[nodiscard]] int getIndex() const;
    [[nodiscard]] int getCharge() const;
    [[nodiscard]] int getLength() const;
    [[nodiscard]] int getNoCal() const;
    [[nodiscard]] float getMz() const;
    [[nodiscard]] float getIRT() const;
    [[nodiscard]] float getSRT() const;
    [[nodiscard]] float getLibQValue() const;

    [[nodiscard]] std::vector<Product> getFragments() const;

private:

    int m_index;
    int m_charge;
    int m_length;
    int m_noCal;
    float m_mz;
    float m_iRT;
    float m_sRT;
    float m_libQvalue;

    std::vector<Product> m_fragments;

};


#endif //PYTHIADIACPP_PEPTIDE_H
