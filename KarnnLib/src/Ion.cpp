//
// Created by anichols on 8/23/23.
//

#include "Ion.h"

Ion::Ion(
        int _type,
        int _index,
        int _loss,
        int _charge,
        double _mz
)
: type(_type)
, index(_index)
, loss(_loss)
, charge(_charge)
, mz(_mz)
{}