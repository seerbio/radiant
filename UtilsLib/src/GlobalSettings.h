//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include "UtilsLib_Exports.h"

#include <QPointF>
#include <QString>
#include <QVector>

using ProteinId = int;
using ProteinSequence = QString;
using ScanNumber = int;
using ScanPoint = QPointF;
using ScanPoints = QVector<QPointF>;

class UTILSLIB_EXPORTS GlobalSettings {

public:

    const QString MODIFICATION_STRING_FORMAT = QStringLiteral("%1%2%3");
    const QChar MODIFICATION_INTERNAL_SEP = '|';
    const QChar SEPARATOR = ';';
    const int HASHING_PRECISION = 3;
    const double ISO_DIFF = 1.0031;
    const int MAX_CHARGE_TANDEM_DEISOTOPING = 2;
    const int POLYNOMIAL_DEGREE = 11;
    const int ROUNDING_PRECISION = 12;
    const double STDEV_MULTIPLIER = 2.5;

    static QString VERSION();

};


const extern UTILSLIB_EXPORTS GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
