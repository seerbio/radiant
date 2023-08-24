//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_QUANTENTRY_H
#define PYTHIADIACPP_QUANTENTRY_H

#include "LibraryCommon.h"

#include "cpu_info.h"


class Fragment {

public:
    mutable int index; // index of the fragment in the library (i.e. its number among the fragments of the precursor)
    mutable float quantity[QVAL::qN];
    mutable float corr;
};

class PrecursorEntry {
public:
    int index; // index of the precursor in the spectral library
    int run_index;
    mutable bool decoy_found;
    mutable int apex, peak, best_fragment, peak_width;
    mutable float RT, RT_start, RT_stop, iRT, predicted_RT, predicted_iRT, one_over_K0, best_fr_mz, profile_qvalue, qvalue, lib_qvalue, dratio, decoy_qvalue, protein_qvalue, pg_qvalue, gg_qvalue, quantity, ratio, level, norm, quality;
    mutable float pg_quantity, pg_norm, gene_quantity, gene_norm, max_lfq, max_lfq_unique;
    mutable float evidence, decoy_evidence, ms1_corr, ms1_area, cscore, decoy_cscore;
    friend inline bool operator < (const PrecursorEntry &left, const PrecursorEntry &right) { return left.index < right.index; }
};

class QuantEntry {
public:
    int index = -1, window, fr_n = 0; // index must be initialised with -1
    mutable PrecursorEntry pr;
    mutable Fragment fr[TopF];

    QuantEntry() { memset(&pr, 0, sizeof(PrecursorEntry)); memset(fr, 0, TopF * sizeof(Fragment)); }

    template <class F> void write(F &out) { out.write((char*)this, sizeof(QuantEntry)); }

    template <class F> void read(F &in) { in.read((char*)this, sizeof(QuantEntry)); }

    friend inline bool operator < (const QuantEntry &left, const QuantEntry &right) { return left.index < right.index; }
};

class DecoyEntry {
public:
    int index = -1;
    float qvalue = 1.0, lib_qvalue = 0.0, dratio = 1.0, cscore = 0.0; // constructor below must have these intialised explicitly

    DecoyEntry() {  }
    DecoyEntry(int _index, float _qvalue, float _lib_qvalue, float _dratio, float _cscore) { index = _index; qvalue = _qvalue; lib_qvalue = _lib_qvalue; dratio = _dratio; cscore = _cscore; }

    template <class F> void write(F &out) {
        out.write((char*)this, sizeof(DecoyEntry));
    }

    template <class F> void read(F &in) {
        in.read((char*)this, sizeof(DecoyEntry));
    }

    friend inline bool operator < (const DecoyEntry &left, const DecoyEntry &right) { return left.index < right.index; }
};


#endif //PYTHIADIACPP_QUANTENTRY_H
