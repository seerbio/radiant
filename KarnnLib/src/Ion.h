//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_ION_H
#define PYTHIADIACPP_ION_H

#include "KarnnLib_Exports.h"


class KARNNLIB_EXPORTS Ion {

public:
    char type;
    char index;
    char loss;
    char charge; // index - number of AAs to the left from the cleavage site
    float mz;

    Ion() = default;

    Ion(
        int _type,
        int _index,
        int _loss,
        int _charge,
        double _mz
        );

    template<class T>
    void init(T &other) {
        type = other.type & 3;
        index = other.index;
        loss = other.loss;
        charge = other.charge;
        mz = other.mz;
    }

    bool operator == (const Ion &other) {
        return index == other.index && type == other.type && charge == other.charge && loss == other.loss;
    }
};


#endif //PYTHIADIACPP_ION_H
