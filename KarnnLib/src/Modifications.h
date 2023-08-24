//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_MODIFICATIONS_H
#define PYTHIADIACPP_MODIFICATIONS_H

#include "KarnnLib_Exports.h"

#include <string>
#include <vector>


struct KARNNLIB_EXPORTS MOD {

    std::string name;
    std::string aas;
    float mass = 0.0;
    int label = 0;

    MOD() = default;

    MOD(
        const std::string &_name,
        float _mass,
        int _label = 0
    )
    : name(_name)
    , mass(_mass)
    , label(_label) {}

    void init(
            const std::string &_name,
            const std::string &_aas,
            float _mass,
            int _label = 0
    ) {
        name = _name;
        aas = _aas;
        mass = _mass;
        label = _label;
    }

    friend bool operator < (const MOD &left, const MOD &right) { return left.name < right.name || (left.name == right.name && left.mass < right.mass); }
};

namespace ModificationsNamespace {

    std::vector<MOD> FixedMods, VarMods;

    extern const KARNNLIB_EXPORTS std::vector<MOD> Modifications;


}



#endif //PYTHIADIACPP_MODIFICATIONS_H
