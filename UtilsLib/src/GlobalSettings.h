//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include "UtilsLib_Exports.h"

#include <QPointF>
#include <QString>
#include <QVector>

using FragLibIonPeptideId = int;
using IsolationWindowKey = QString;
using ModificationMass = double;
using MsLevel = int;
using MZION = double;
using MzMin = double;
using MzMax = double;
using MzTargetKey = QString;
using PeptideId = int;
using PeptideString = QString;
using ProteinId = int;
using ProteinSequence = QString;
using ResidueIndex = int;
using ScanNumber = int;
using ScanPoint = QPointF;
using ScanPoints = QVector<QPointF>;
using TandemScansIndex = int;
using TARGETMZ = double;
using UniqueMsInfoScanKey = QString;
using UniqueHashedMzAndTarget = QString;

class UTILSLIB_EXPORTS GlobalSettings {

public:

    const QChar COMMA = ',';
    const QString DOT_CACHE = QStringLiteral(".cache");
    const QString DOT_FRAGLIB = QStringLiteral(".fragLib");
    const QString DOT_PEPLIB = QStringLiteral(".pepLib");
    const int HASHING_PRECISION = 3;
    const double ISO_DIFF = 1.0031;
    const int MAX_CHARGE_TANDEM_DEISOTOPING = 2;
    const QChar MODIFICATION_INTERNAL_SEP = '|';
    const QString MODIFICATION_STRING_FORMAT = QStringLiteral("%1%2%3");
    const QString NEWLINE = "\n";
    const QString NONE = QStringLiteral("NONE");
    const int POLYNOMIAL_DEGREE = 11;
    const int ROUNDING_PRECISION = 12;
    const QChar SEPARATOR = ';';
    const double STDEV_MULTIPLIER = 2.5;
    const QString TAB = "\t";

    static QString VERSION();
};


const extern UTILSLIB_EXPORTS GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
