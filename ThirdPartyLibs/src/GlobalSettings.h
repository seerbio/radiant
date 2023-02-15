//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include <string>
#include <vector>


class  GlobalSettings {

public:

    const std::string MODIFICATION_STRING_FORMAT = "%1%2%3";
    const char MODIFICATION_INTERNAL_SEP = '|';
    const char SEPARATOR = ';';
    const char COMMA = ',';
    const std::string TAB = "\t";
    const std::string NEWLINE = "\n";
    const int HASHING_PRECISION = 3;
    const double ISO_DIFF = 1.0031;
    const int MAX_CHARGE_TANDEM_DEISOTOPING = 2;
    const int POLYNOMIAL_DEGREE = 11;
    const int ROUNDING_PRECISION = 12;
    const double STDEV_MULTIPLIER = 2.5;

    static std::string VERSION();

};


const extern GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
