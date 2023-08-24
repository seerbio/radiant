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


namespace LibraryCommon {

    static double AA[256];

    extern const int KARNNLIB_EXPORTS fTypeB;
    extern const int KARNNLIB_EXPORTS fTypeY;
    extern const int KARNNLIB_EXPORTS fExclude;

    extern const int KARNNLIB_EXPORTS MaxF;
    extern const bool KARNNLIB_EXPORTS ReverseDecoys;
    extern const int KARNNLIB_EXPORTS MaxRecCharge;
    extern const int KARNNLIB_EXPORTS MaxRecLoss;
    extern const double KARNNLIB_EXPORTS GlobalMassAccuracy;
    extern const double KARNNLIB_EXPORTS GlobalMassAccuracyMs1;
    extern const double KARNNLIB_EXPORTS GeneratorAccuracy;

    extern const double KARNNLIB_EXPORTS PROTON;
    extern const double KARNNLIB_EXPORTS OH;
    extern const double KARNNLIB_EXPORTS C13delta;

    extern const std::vector<double> KARNNLIB_EXPORTS Loss;


    template <class T>
    double KARNNLIB_EXPORTS sum(T *x, int n) {
        double ex = 0.0;
        for (int i = 0; i < n; i++) ex += x[i];
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

    int closing_bracket(
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

    int KARNNLIB_EXPORTS peptide_length(const std::string &name);

    std::string KARNNLIB_EXPORTS get_aas(const std::string &name);

    std::string KARNNLIB_EXPORTS getLoss(int lossIndex);

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

};


#endif //PYTHIADIACPP_LIBRARYCOMMON_H
