//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_PEPTIDE_H
#define PYTHIADIACPP_PEPTIDE_H

#include "KarnnLib_Exports.h"

#include "Product.h"
#include "ReadWriteCommonFunctons.h"

#include <vector>

class KARNNLIB_EXPORTS Peptide {

public:

    int index = 0;
    int charge = 0;
    int length = 0;
    int no_cal = 0;
    float mz = 0.0;
    float iRT = 0.0;
    float sRT = 0.0;
    float lib_qvalue = 0.0;

    std::vector<Product> fragments;

    void init(float _mz, float _iRT, int _charge, int _index) {
        mz = _mz;
        iRT = _iRT;
        sRT = 0.0;
        charge = _charge;
        index = _index;
    }

    inline void free() {
        std::vector<Product>().swap(fragments);
    }

    template <class F>
    void write(F &out) {
        out.write((char*)&index, sizeof(int));
        out.write((char*)&charge, sizeof(int));
        out.write((char*)&length, sizeof(int));

        out.write((char*)&mz, sizeof(float));
        out.write((char*)&iRT, sizeof(float));
        out.write((char*)&sRT, sizeof(float));

        ReadWriteCommonFunctons::write_vector(out, fragments);
    }

    template <class F>
    void read(F &in) {
        in.read((char*)&index, sizeof(int));
        in.read((char*)&charge, sizeof(int));
        in.read((char*)&length, sizeof(int));

        in.read((char*)&mz, sizeof(float));
        in.read((char*)&iRT, sizeof(float));
        in.read((char*)&sRT, sizeof(float));

        ReadWriteCommonFunctons::read_vector(in, fragments);
    }

};


#endif //PYTHIADIACPP_PEPTIDE_H
