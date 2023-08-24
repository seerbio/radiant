//
// Created by anichols on 8/23/23.
//

#include "Product.h"


Product::Product(
        float _mz,
        float _height,
        int _charge)
: mz(_mz)
, height(_height)
, charge(_charge)
{}

Product::Product(
        float _mz,
        float _height,
        int _charge,
        int _type,
        int _index,
        int _loss
        )
: mz(_mz)
, height(_height)
, charge(_charge)
, type(_type)
, index(_index)
, loss(_loss)
{}

void Product::init(
        float _mz,
        float _height,
        int _charge
) {
    mz = _mz;
    height = _height;
    charge = _charge;
}

int Product::ion_code() {
    return (((((int)type) * 20 + (int)charge)) * (loss_other + 1) + (int)loss) * 100 + (int)index + 1;
}