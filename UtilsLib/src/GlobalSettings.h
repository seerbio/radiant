//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include "UtilsLib_Exports.h"

#include "ErrorUtils.h"
#include "PointDF.h"
#include "PointFF.h"

#include <QPointF>
#include <QString>
#include <QVector>

using namespace Error;

using Charge = int;
using Coors = QVector<double>;
using CosineSimSum = double;
using DecoyRatio = double;
using DiscScore = double;
using FilePath = QString;
using FragLibIonPeptideId = int;
using FrameIndex = int;
using FrequencyPercent = double;
using Id = int;
using Index = int;
using Intensity = double;
using IonIndex = int;
using IonLabel = QString;
using IonLabels = QStringList;
using IonMobilityIndex = int;
using IonType = QString;
using IRT = float;
using IsolationWindowKey = QString;
using Mass = double;
using ModificationMass = double;
using MonoOffset = int;
using MsLevel = int;
using MzHashed = int;
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
using QValue = double;
using ReScore = double;
using ResidueIndex = int;
using ScanNumber = int;
using ScanNumberIndex = int;
using ScanPoint = PointFF;
using ScanPoints = QVector<ScanPoint>;
using ScanTime = float;
using Score = double;
using TandemScansIndex = int;
using TARGETMZ = double;
using Tranche = int;
using MzTargetKey = QString;
using UniqueHashedMzAndTarget = QString;
using Value = double;
using XVal = double;
using YVal = double;


class UTILSLIB_EXPORTS GlobalSettings {

public:

    const QString AAs = "GAVLIFMPWSCTYHKRQEND";
    const QString MutateAAto = "LLLVVLLLLTSSSSLLNDQE";

    const QChar COMMA = ',';
    const QString DOT_CACHE = QStringLiteral(".cache");
    const QString DOT_CAL = QStringLiteral(".cal");
    const QString DOT_CSV = QStringLiteral(".csv");
    const QString DOT_FASTA = QStringLiteral(".fasta");
    const QString DOT_FRAGLIB_DF = QStringLiteral(".fragLibDF");
    const QString DOT_FRAGLIB_FF = QStringLiteral(".fragLibFF");
    const QString DOT_LIB = QStringLiteral(".lib");
    const QString DOT_MAT = QStringLiteral(".mat");
    const QString DOT_PRQ_DF = QStringLiteral(".prqDF");
    const QString DOT_PEPLIB = QStringLiteral(".pepLib");
    const QString DOT_PSM = QStringLiteral(".psm");
    const QString FAILED_SHUFFLE = QStringLiteral("Failed Shuffle");
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
    const QString PRQ_FILE_EXTENSION = QStringLiteral("prqFF");
    const QString PSM_FILE_EXTENSION = QStringLiteral("psm");
    const QString PSM_SCORED_FILE_EXTENSION = QStringLiteral("scored");
    const QString PYTHIA_FILE_EXTENSION = QStringLiteral("pythia");
    const QString PYTHIA_CAL_FILE_EXTENSION = QStringLiteral("pythiaCAL");
    const QString PYTHIA_DIA_FILE_EXTENSION = QStringLiteral("pythiaDIA");

    const QString DOT_CACHED_FILE_EXTENSION = QStringLiteral(".cached");
    const QString DOT_CSV_FILE_EXTENSION = QStringLiteral(".csv");
    const QString DOT_FASTA_FILE_EXTENSION = QStringLiteral(".fasta");
    const QString DOT_HDF_FILE_EXTENSION = QStringLiteral(".hdf");
    const QString DOT_MZML_FILE_EXTENSION = QStringLiteral(".mzml");
    const QString DOT_PRQ_FILE_EXTENSION = QStringLiteral(".prq");
    const QString DOT_PSM_FILE_EXTENSION = QStringLiteral(".psm");
    const QString DOT_PSM_SCORED_FILE_EXTENSION = QStringLiteral(".scored");
    const QString DOT_PYTHIA_FILE_EXTENSION = QStringLiteral(".pythia");
    const QString DOT_PYTHIA_CAL_FILE_EXTENSION = QStringLiteral(".pythiaCAL");
    const QString DOT_PYTHIA_DIA_FILE_EXTENSION = QStringLiteral(".pythiaDIA");

    const QString Y_IONS = QStringLiteral("yIons");
    const QString B_IONS = QStringLiteral("bIons");
    const QString Y2_IONS = QStringLiteral("y2Ions");
    const QString B2_IONS = QStringLiteral("b2Ions");
    const QString A_IONS = QStringLiteral("aIons");
    const QString Y_NH3_IONS = QStringLiteral("yNH3");
    const QString Y_H2O_IONS = QStringLiteral("yH20");
    const QString B_NH3_IONS= QStringLiteral("bNH3");
    const QString B_H2O_IONS = QStringLiteral("bH20");
    const QString PRECURSOR_IONS = QStringLiteral("precursorIons");

    const double TIGHT_1_FRACTION = 0.45;
    const double TIGHT_2_FRACTION = 0.2;

    const int NUMBER_OF_THE_BEAST = 666;

    const QString MS1Key = QStringLiteral("MS1Key");

    static QString VERSION();
};


const extern UTILSLIB_EXPORTS GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
