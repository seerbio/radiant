//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include "UtilsLib_Exports.h"

#include <QPointF>
#include <QString>
#include <QVector>

using Charge = int;
using Coors = QVector<double>;
using DiscScore = double;
using FragLibIonPeptideId = int;
using FrameIndex = int;
using Id = int;
using Index = int;
using Intensity = double;
using IonLabel = QString;
using IonLabels = QStringList;
using IonMobilityIndex = int;
using IsolationWindowKey = QString;
using Mass = double;
using ModificationMass = double;
using MonoOffset = int;
using MS2Ion = QPointF;
using MsLevel = int;
using MZION = double;
using MzMin = double;
using MzMax = double;
using MzTargetKey = QString;
using NominalMzMass = int;
using Occurrence = int;
using PeakIntegrationIndexes = QPair<int, int>;
using PeptideId = int;
using PeptideSequenceChargeKey = QString;
using PeptideString = QString;
using PeptideStringWithMods = QString;
using DiffPPM = double;
using ProteinId = int;
using ProteinSequence = QString;
using ReScore = double;
using ResidueIndex = int;
using ScanNumber = int;
using ScanNumberIndex = int;
using ScanPoint = QPointF;
using ScanPoints = QVector<QPointF>;
using ScanTime = double;
using Score = double;
using TandemScansIndex = int;
using TARGETMZ = double;
using Tranche = int;
using UniqueMsInfoScanKey = QString;
using UniqueHashedMzAndTarget = QString;
using Value = double;
using XICPoint = QPair<ScanNumber, Intensity>;
using XICPoints = QMap<ScanNumber, Intensity>;

class UTILSLIB_EXPORTS GlobalSettings {

public:

    const QChar COMMA = ',';
    const QString DOT_CACHE = QStringLiteral(".cache");
    const QString DOT_CAL = QStringLiteral(".cal");
    const QString DOT_CSV = QStringLiteral(".csv");
    const QString DOT_FASTA = QStringLiteral(".fasta");
    const QString DOT_FRAGLIB = QStringLiteral(".fragLib");
    const QString DOT_LIB = QStringLiteral(".lib");
    const QString DOT_MAT = QStringLiteral(".mat");
    const QString DOT_PRQ = QStringLiteral(".prq");
    const QString DOT_PEPLIB = QStringLiteral(".pepLib");
    const QString DOT_PSM = QStringLiteral(".psm");
    const int HASHING_PRECISION = 3;
    const double ISO_DIFF = 1.003355;
    const int MAX_CHARGE_TANDEM_DEISOTOPING = 2;
    const QChar MODIFICATION_INTERNAL_SEP = '|';
    const QString MODIFICATION_STRING_FORMAT = QStringLiteral("%1%2%3");
    const QString NEWLINE = "\n";
    const QString NONE = QStringLiteral("NONE");
    const int POLYNOMIAL_DEGREE = 11;
    const int ROUNDING_PRECISION = 12;
    const QChar SEPARATOR = ';';
    const double STDEV_MULTIPLIER = 3.0;
    const QString TAB = "\t";

    const QString CACHED_FILE_EXTENSION = QStringLiteral("cached");
    const QString CSV_FILE_EXTENSION = QStringLiteral("csv");
    const QString FASTA_FILE_EXTENSION = QStringLiteral("fasta");
    const QString HDF_FILE_EXTENSION = QStringLiteral("hdf");
    const QString MZML_FILE_EXTENSION = QStringLiteral("mzml");
    const QString PRQ_FILE_EXTENSION = QStringLiteral("prq");
    const QString PSM_FILE_EXTENSION = QStringLiteral("psm");
    const QString PSM_SCORED_FILE_EXTENSION = QStringLiteral("scored");
    const QString PYTHIA_FILE_EXTENSION = QStringLiteral("pythia");

    static QString VERSION();
};


const extern UTILSLIB_EXPORTS GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
