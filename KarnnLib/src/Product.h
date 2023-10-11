//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_PRODUCT_H
#define PYTHIADIACPP_PRODUCT_H

#include "KarnnLib_Exports.h"

#include <string>

enum {
    loss_none,
    loss_H2O,
    loss_NH3,
    loss_CO,
    loss_N,
    loss_other
};


class KARNNLIB_EXPORTS Product {

public:

    float mz = 0.0;
    float height = 0.0;
    char charge = 0;
    char type = 0;
    char index = 0;
    char loss = 0;

    Product() = default;

    Product(
            float _mz,
            float _height,
            int _charge
            );

    Product(
            float _mz,
            float _height,
            int _charge,
            int _type,
            int _index,
            int _loss
            );

    static std::string getLoss(int lossIndex) {

        switch (lossIndex) {

            case 0:
                return "";
            case 1:
                return "-H2O";
            case 2:
                return "-NH3";
            case 3:
                return "a";
            case 4:
                return "N";
            case 5:
                return "other";
            default:
                return "";
        }
    }

    friend inline bool operator < (const Product &left, const Product &right) {return left.mz < right.mz;}

    void init(
            float _mz,
            float _height,
            int _charge
            );

    int ion_code();

};


#endif //PYTHIADIACPP_PRODUCT_H
