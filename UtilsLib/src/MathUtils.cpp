//
// Created by anichols on 11/07/2021.
//

#include "MathUtils.h"




namespace MathUtilsConstants {

    const double PI = 3.14159265359;


}//NAMESPACE


double MathUtils::calculatePPM(double val, double ppmTolerance) {
    return (val * ppmTolerance) / 1e6;
}


double MathUtils::factorial(int n) {

    if (n < 0 ) {
        return 0;
    }

    const double f = !n ? static_cast<double>(1) : n * factorial(n - 1);

    return f < 0 ? std::numeric_limits<double>::max() : f;
}
