//
// Created by anichols on 8/23/23.
//

#include "LibraryCommon.h"

#include "cpu_info.h"
#include "Modifications.h"
#include "Product.h"

#include <QString>
#include <QDebug>


namespace LibraryCommon {

    const double PROTON = 1.007825035;
    const double OH = 17.003288;
    const double C13delta = 1.003355;

    const int fTypeB = 1 << 0;
    const int fTypeY = 1 << 1;
    const int fExclude = 1 << 6;

    const std::vector<double> KARNNLIB_EXPORTS Loss = { 0.0, 18.011113035, 17.026548, 27.994915, 0.0, 0.0 };

    const int MaxF = std::numeric_limits<int>::max();
    const bool ReverseDecoys = false;
    const int MaxRecCharge = 19;
    const int MaxRecLoss = loss_N;
    const double GlobalMassAccuracy = 20.0 / 1000000.0;
    const double GlobalMassAccuracyMs1 = 20.0 / 1000000.0;
    const double GeneratorAccuracy = 5.0 / 1000000.0;

    int closing_bracket(
            const std::string &name,
            char symbol,
            int pos
            ) {

        char close = (symbol == '(' ? ')' : ']');

        int end;
        int par;

        for (end = pos + 1, par = 1; end < name.size(); end++) {

            char s = name[end];

            if (s == close) {
                par--;
                if (!par) {
                    break;
                }
            }

            else {
                if (s == symbol) {
                    par++;
                }
            }
        }

        return end;
    }

    int peptide_length(const std::string &name) {

        int n = 0;

        for (int i = 0; i < name.size(); i++) {

            char symbol = name[i];

            if (symbol < 'A' || symbol > 'Z') {

                if (symbol != '(' && symbol != '[') {
                    continue;
                }

                i = closing_bracket(name, symbol, i);

            }

            n++;
        }
        return n;
    }

    std::string get_aas(const std::string &name) {
        int i, end, par;
        std::string result;

        for (i = 0; i < name.size(); i++) {
            char symbol = name[i];
            if (symbol < 'A' || symbol > 'Z') {
                if (symbol != '(' && symbol != '[') continue;
                i = closing_bracket(name, symbol, i);
                continue;
            }
            result.push_back(symbol);
        }
        return result;
    }

    std::string getLoss(int lossIndex) {

        switch (lossIndex) {

            case loss_none:
                return "";
            case loss_H2O:
                return "-H2O";
            case loss_NH3:
                return "-NH3";
            case loss_CO:
                return "a";
            case loss_N:
                return "N";
            case loss_other:
                return "other";
            default:
                return "";
        }

    }

    std::vector<Ion> KARNNLIB_EXPORTS recognise_fragments(
            const std::vector<double> &sequence,
            const std::vector<Product> &fragments,
            float &gen_acc,
            int pr_charge,
            bool full_spectrum,
            int max_charge,
            int loss_cap
                    ) {

        double max_gen_acc = 0.0;
        Lock gen_acc_lock;

        int i, j, cnt, tot, charge, loss, index;
        std::vector<Ion> result(fragments.size());
        std::vector<Ion> v;
        double delta, min, s = sum(&(sequence[0]), sequence.size());
        gen_acc = GeneratorAccuracy;

        anew:
        cnt = tot = index = 0, charge = 1;
        loss = loss_none;
        for (i = 0; i < result.size(); i++) result[i].charge = 0;

        start:
        v = generate_fragments(sequence, charge, loss, &cnt);
        if (charge == 1) {
            v.push_back(Ion(fTypeY, sequence.size(), loss, pr_charge, (s + OH + (1.0 + double(pr_charge)) * PROTON - Loss[loss]) / static_cast<double>(pr_charge))); // non-fragmented
        }

        for (i = 0; i < fragments.size(); i++) {
            if (result[i].charge) continue;
            double margin = fragments[i].mz * gen_acc;
            for (j = 0, min = margin; j < v.size(); j++) {
                delta = std::abs(v[j].mz - fragments[i].mz);
                if (delta < min) min = delta, index = j;
            }
            if (min < margin) {
                tot++;
                result[i].init(v[index]);
                if (tot >= fragments.size()) goto stop;
            }
        }

        loss++;
        if (loss == loss_cap) loss = loss_none, charge++;
        if (charge <= max_charge && tot < fragments.size()) goto start;

        stop:
        if (tot < fragments.size()) {
            if (full_spectrum) return result;
            gen_acc *= 2.0;

            if (gen_acc > max_gen_acc) {
                while (!gen_acc_lock.set()) {}
                max_gen_acc = gen_acc;
                gen_acc_lock.free();
            }
            v.clear();
            goto anew;
        }

        return result;
    }

    std::vector<Ion> KARNNLIB_EXPORTS recognise_fragments(
            const std::string &name,
            const std::vector<Product> &fragments,
            float &gen_acc,
            int pr_charge,
            bool full_spectrum,
            int max_charge,
            int loss_cap
            ) {

        int i;

        if (fragments.size()) if (fragments[0].type & 3) {

            std::vector<Ion> result(fragments.size());

            for (i = 0; i < result.size(); i++) {
                if (fragments[i].type & 3) result[i].init(fragments[i]); else break;
            }
            if (i == result.size()) {
                return result;
            }
        }

        auto sequence = get_sequence(name);
        return recognise_fragments(sequence, fragments, gen_acc, pr_charge, full_spectrum, max_charge, loss_cap);
    }


    std::vector<double> get_sequence(
            const std::string &name,
            int *no_cal
    ) {
        int i, j;
        double add = 0.0;
        std::vector<double> result;
        if (no_cal != NULL) *no_cal = 0;

        for (i = 0; i < name.size(); i++) {

            char symbol = name[i];

            if (symbol < 'A' || symbol > 'Z') {

                if (symbol != '(' && symbol != '[') {
                    continue;
                }

                i++;

                int end = closing_bracket(name, symbol, i);

                if (end >= name.size()) {
                    qDebug() << "incorrect peptide name format: ";
                }

                std::string mod = name.substr(i, end - i);
                for (j = 0; j < ModificationsNamespace::Modifications.size(); j++) {

                    if (ModificationsNamespace::Modifications[j].name == mod) {

                        if (!result.size()) {
                            add += ModificationsNamespace::Modifications[j].mass;
                        }

                        else {
                            result[result.size() - 1] += ModificationsNamespace::Modifications[j].mass;
                        }

                        if (
                            no_cal != nullptr
                            && ModificationsNamespace::Modifications[j].mass < 4.5
                            && ModificationsNamespace::Modifications[j].mass > 0.0
                            && std::abs(ModificationsNamespace::Modifications[j].mass - static_cast<double>(floor(ModificationsNamespace::Modifications[j].mass + 0.5)) < 0.1)
                            ) {
                            *no_cal = 1;
                        }

                        break;
                    }
                }

                if (j == ModificationsNamespace::Modifications.size()){
                    qDebug() << QStringLiteral("unknown modification: ");
                }

                i = end;
                continue;
            }

            result.push_back(AA[symbol] + add);
            add = 0.0;
        }
        return result;
    }

    std::vector<Ion>  KARNNLIB_EXPORTS generate_fragments(
            const std::vector<double> &sequence,
            int charge,
            int loss,
            int *cnt,
            int min_aas
            ) {

        int i;
        double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());

        std::vector<Ion> result;

        for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
            curr += sequence[i];
            double b = curr;
            double y = s - curr + PROTON + OH;

            if (i + 1 >= min_aas) {
                result.emplace_back(fTypeB, i + 1, loss, charge, (b + c * PROTON - Loss[loss]) / c), (*cnt)++;
            }

            if (sequence.size() - i - 1 >= min_aas) {
                result.emplace_back(fTypeY, i + 1, loss, charge, (y + c * PROTON - Loss[loss]) / c), (*cnt)++;
            }
        }
        return result;

    }


}