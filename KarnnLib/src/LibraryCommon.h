//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_LIBRARYCOMMON_H
#define PYTHIADIACPP_LIBRARYCOMMON_H

#include "KarnnLib_Exports.h"

#include "Ion.h"
#include "Product.h"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

const int TopF = 6;
const int auxF = 12;
const int nnS = 2;
const int nnW = (2 * nnS) + 1;
const int qL = 3;
const int QSL[qL] = { 1, 3, 5 };

enum QVAL {
    qTotal,
    qFiltered,
    qN
};

enum LIB {
    libPr, libCharge, libPrMz,
    libiRT, libFrMz, libFrI,
    libPID, libPN, libGenes, libPT, libIsDecoy,
    libFrCharge, libFrType, libFrNumber, libFrLoss,
    libQ, libEG, libFrExc, libCols
};

enum FEATURES {
    pTimeCorr,
    pLocCorr, pMinCorr, pTotal, pCos, pCosCube, pMs1TimeCorr, pNFCorr, pdRT, pResCorr, pResCorrNorm, pTightCorrOne, pTightCorrTwo, pShadow, pHeavy,
    pMs1TightOne, pMs1TightTwo,
    pMs1Iso, pMs1IsoOne, pMs1IsoTwo,
    pMs1Ratio,
    pBestCorrDelta, pTotCorrSum,
    pMz, pCharge, pLength, pFrNum,
    pMods,
    pAAs,
    pRT = pAAs + 20,
    pAcc,
    pRef = pAcc + TopF,
    pSig = pRef + auxF - 1,
    pCorr = pSig + TopF,
    pShadowCorr = pCorr + auxF,
    pShape = pShadowCorr + TopF,
    pQLeft = pShape + nnW,
    pQRight = pQLeft + qL,
    pQPos = pQRight + qL,
    pQNFCorr = pQPos + qL,
    pQCorr = pQNFCorr + qL,
    pN = pQCorr + qL
};


namespace LibraryCommon {

    static double AA[256];
    static int AA_index[256];
    static std::vector<std::pair<std::string, std::string> > UniMod;
    static std::vector<std::pair<std::string, int> > UniModIndex;
    static std::vector<int> UniModIndices;

    extern const KARNNLIB_EXPORTS char * AAs;
    extern const KARNNLIB_EXPORTS char * MutateAAto;

    extern const int KARNNLIB_EXPORTS fTypeB;
    extern const int KARNNLIB_EXPORTS fTypeY;
    extern const int KARNNLIB_EXPORTS fExclude;

    extern const int KARNNLIB_EXPORTS MaxF;
    extern const bool KARNNLIB_EXPORTS ReverseDecoys;
    extern const int KARNNLIB_EXPORTS MaxRecCharge;
    extern const int KARNNLIB_EXPORTS MaxRecLoss;
    extern const bool KARNNLIB_EXPORTS ForceFragRec;
    extern const bool KARNNLIB_EXPORTS ExportRecFragOnly;
    extern const double KARNNLIB_EXPORTS GlobalMassAccuracy;
    extern const double KARNNLIB_EXPORTS GlobalMassAccuracyMs1;
    extern const double KARNNLIB_EXPORTS GeneratorAccuracy;
    extern const double KARNNLIB_EXPORTS SptxtMassAccuracy;

    extern const double KARNNLIB_EXPORTS PROTON;
    extern const double KARNNLIB_EXPORTS OH;
    extern const double KARNNLIB_EXPORTS C13delta;

    extern const std::vector<double> KARNNLIB_EXPORTS Loss;

    template <class T>
    double KARNNLIB_EXPORTS sum(T *x, int n) {
        double ex = 0.0;
        for (int i = 0; i < n; i++) {
            ex += x[i];
        }
        return ex;
    }

    template <class T>
    double KARNNLIB_EXPORTS sum(std::vector<T> &x) {
        return sum(&(x[0]), x.size());
    }

    template <class T>
    double KARNNLIB_EXPORTS mean(T * x, int n) {
        assert(n >= 1);
        return sum(x, n) / (double)n;
    }

    template <class T>
    double KARNNLIB_EXPORTS mean(std::vector<T> &x) {
        int n = x.size();
        assert(n >= 1);
        return sum(x) / (double)n;
    }

    template <class T>
    double KARNNLIB_EXPORTS var(T * x, int n) {
        double ex = mean(x, n), s2 = 0.0;
        for (int i = 0; i < n; i++) s2 += Sqr(x[i] - ex);
        if (n > 1) s2 /= (double)(n - 1);
        return s2;
    }

    template <class T>
    double KARNNLIB_EXPORTS var(std::vector<T> &x) {
        return var(&(x[0]), x.size());
    }

    template <class Tx, class Ty>
    double KARNNLIB_EXPORTS cos(Tx * x, Ty * y, int n) {

        const double E = std::numeric_limits<double>::min();

        double s2x = 0.0, s2y = 0.0, r = 0.0;
        for (int i = 0; i < n; i++) s2x += Sqr(x[i]), s2y += Sqr(y[i]), r += x[i] * y[i];
        if (s2x < E || s2y < E) return 0.0;
        return r / sqrt(s2x * s2y);
    }

    template <class Tx, class Ty>
    double KARNNLIB_EXPORTS scalar(Tx * x, Ty * y, int n) {
        double r = 0.0;
        for (int i = 0; i < n; i++) r += x[i] * y[i];
        return r;
    }

    template <class Tx, class Ty>
    double KARNNLIB_EXPORTS scalar(Tx * x, Ty * y, bool * mask, int n) {
        double r = 0.0;
        for (int i = 0; i < n; i++) if (mask[i]) r += x[i] * y[i];
        return r;
    }

    int KARNNLIB_EXPORTS y_delta_composition_index(char aa);

    int KARNNLIB_EXPORTS y_delta_charge_index();

    int KARNNLIB_EXPORTS y_delta_index(char aa, int shift);

    int KARNNLIB_EXPORTS y_cterm_index(char aa, int shift);

    int KARNNLIB_EXPORTS y_nterm_index(int shift);

    int KARNNLIB_EXPORTS y_delta_size();

    double KARNNLIB_EXPORTS y_ratio(
            int i,
            int charge,
            const std::string &aas
            );

    void KARNNLIB_EXPORTS y_scores(
            std::vector<double> &s,
            int charge,
            const std::string &aas
            );

    void KARNNLIB_EXPORTS to_eg(std::string &eg, const std::string &name);

    std::string KARNNLIB_EXPORTS to_eg(const std::string &name);

    void KARNNLIB_EXPORTS to_exp(std::vector<double> &s);

    int KARNNLIB_EXPORTS y_loss_aa(char aa, bool NH3);

    int KARNNLIB_EXPORTS y_loss(bool NH3);

    int KARNNLIB_EXPORTS y_loss_size();

    void KARNNLIB_EXPORTS y_loss_scores(
            std::vector<double> &s,
            const std::string &aas,
            bool NH3
            );

    int KARNNLIB_EXPORTS closing_bracket(
            const std::string &name,
            char symbol,
            int pos
    );

    std::vector<Ion>  KARNNLIB_EXPORTS generate_fragments(
            const std::vector<double> &sequence,
            int charge,
            int loss,
            int *cnt,
            int min_aas = 0
                    );

    std::vector<Ion> KARNNLIB_EXPORTS generate_fragments(
            const std::vector<double> &sequence,
            const std::vector<Ion> &pattern,
            int pr_charge
            );

    int KARNNLIB_EXPORTS peptide_length(const std::string &name);

    std::string KARNNLIB_EXPORTS get_aas(const std::string &name);

    std::string KARNNLIB_EXPORTS getLoss(int lossIndex);

    std::vector<Ion> KARNNLIB_EXPORTS generate_all_fragments(
            const std::vector<double> &sequence,
            int charge,
            int max_loss,
            int *cnt
            );

    std::vector<Ion> KARNNLIB_EXPORTS recognise_fragments(
            const std::vector<double> &sequence,
            const std::vector<Product> &fragments,
            float &gen_acc,
            int pr_charge,
            bool full_spectrum = false,
            int max_charge = 19,
            int loss_cap = loss_N
                    );

    std::vector<Ion> KARNNLIB_EXPORTS recognise_fragments(
            const std::string &name,
            const std::vector<Product> &fragments,
            float &gen_acc,
            int pr_charge,
            bool full_spectrum = false,
            int max_charge = 19,
            int loss_cap = loss_N
                    );

    std::vector<double> KARNNLIB_EXPORTS get_sequence(
            const std::string &name,
            int *no_cal = NULL
                    );

    std::vector<bool> KARNNLIB_EXPORTS get_mod_state(
            const std::string &name,
            int length
            );

    std::string get_extension(const std::string &file);

};


#endif //PYTHIADIACPP_LIBRARYCOMMON_H
