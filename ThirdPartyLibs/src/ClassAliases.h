//
// Created by anichols on 2/7/23.
//

#ifndef SPARKDIA_CLASSALIASES_H
#define SPARKDIA_CLASSALIASES_H

#include <vector>


class ScanPoint {

public:

    double x;
    double y;

    ScanPoint(double x, double y) : x(x), y(y) {}
    ~ScanPoint() = default;

};

using ScanNumber = int;
using ScanPoints = std::vector<ScanPoint>;

#endif //SPARKDIA_CLASSALIASES_H
