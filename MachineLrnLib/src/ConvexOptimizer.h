//
// Created by anichols on 9/9/23.
//

#ifndef PYTHIADIACPP_CONVEXOPTIMIZER_H
#define PYTHIADIACPP_CONVEXOPTIMIZER_H

#include "Error.h"
#include "MachineLrnLib_Exports.h"

using namespace Error;


class MACHINELRNLIB_EXPORTS ConvexOptimizer {


public:

    ConvexOptimizer() = default;
    ~ConvexOptimizer() = default;

    Err test();

};


#endif //PYTHIADIACPP_CONVEXOPTIMIZER_H
