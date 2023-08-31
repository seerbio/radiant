//
// Created by anichols on 11/07/2021.
//

#include "MathUtils.h"

#include "ErrorUtils.h"



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


QPointF MathUtils::closestXValPoint(const QVector<QPointF> &vec, double value) {

    assert(!vec.isEmpty());
    auto it = std::min_element(vec.begin(), vec.end(), [value] (const QPointF &a, const QPointF &b) {
        return std::abs(a.x() - value) <= std::abs(b.x() - value);
    });

    if(it == vec.end()) {
        return {-1, -1};
    }

    return *it;
}
