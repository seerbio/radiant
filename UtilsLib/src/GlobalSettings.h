//
// Created by Drucifer on 3/21/2022.
//

#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H

#include "UtilsLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"

#include <QPointF>
#include <QString>
#include <QVector>

using namespace Error;

using Charge = int;
using Coors = QVector<double>;
using CosineSimSum = double;
using DiscScore = double;
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

struct UTILSLIB_EXPORTS MS2Ion {

    double mz = -1.0;
    double intensity = -1.0;
    IRT iRT = -1.0;
    QString ionLabel;

    MS2Ion() = default;

    MS2Ion(
            double mz,
            double intensity,
            QString ionLabel
    )
            :mz(mz)
            , intensity(intensity)
            , ionLabel(std::move(ionLabel))
    {}

    friend QDebug operator<<(QDebug dbg, const MS2Ion& obj) {
        dbg.nospace() << "MS2Ion(" << obj.mz << ", " << obj.intensity << ")";
        return dbg;
    }

    Err getIonLabelInfo(QPair<IonIndex, IonType> *ionInfo) const {

        ERR_INIT

        const int expectedSplitSize = 2;

        IonIndex ionIndex;

        if (ionLabel.contains("p")) {
             *ionInfo = {0, ionLabel};
             ERR_RETURN
        }
        else if (ionLabel.contains("p-H2O")) {
            *ionInfo = {1, ionLabel};
            ERR_RETURN
        }
        else if (ionLabel.contains("p-NH3")) {
            *ionInfo = {2, ionLabel};
            ERR_RETURN
        }
        if (ionLabel.contains("p^2")) {
            *ionInfo = {3, ionLabel};
            ERR_RETURN
        }
        else if (ionLabel.contains("p-H2O^2")) {
            *ionInfo = {4, ionLabel};
            ERR_RETURN
        }
        else if (ionLabel.contains("p-NH3^2")) {
            *ionInfo = {5, ionLabel};
            ERR_RETURN
        }
        else if (ionLabel.contains("y")) {

            QStringList ionLabelSplit;
            if (ionLabel.contains("-")) {

                ionLabelSplit = ionLabel.split("-", QString::SkipEmptyParts);
                e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

                const QString labelFront = ionLabelSplit.front().replace("y", "");

                const QString labelBack = ionLabelSplit.back();
                e = ErrorUtils::toInt(labelFront, &ionIndex); ree

                *ionInfo = {ionIndex, "y-" + labelBack};
                ERR_RETURN
            }
            else if (ionLabel.contains("^")) {
                ionLabelSplit = ionLabel.split("^", QString::SkipEmptyParts);
                e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

                const QString labelFront = ionLabelSplit.front().replace("y", "");

                const QString labelBack = ionLabelSplit.back();
                e = ErrorUtils::toInt(labelFront, &ionIndex); ree

                *ionInfo = {ionIndex, "y^" + labelBack};
                ERR_RETURN
            }
            else {

                QString ionLabelCopy = ionLabel;
                ionLabelCopy = ionLabelCopy.replace("y", "");
                e = ErrorUtils::toInt(ionLabelCopy, &ionIndex); ree
                *ionInfo = {ionIndex, "y"};
                ERR_RETURN
            }
        }

        else if (ionLabel.contains("b")) {

            QStringList ionLabelSplit;
            if (ionLabel.contains("-")) {

                ionLabelSplit = ionLabel.split("-", QString::SkipEmptyParts);
                e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

                const QString labelFront = ionLabelSplit.front().replace("b", "");

                const QString labelBack = ionLabelSplit.back();
                e = ErrorUtils::toInt(labelFront, &ionIndex); ree

                *ionInfo = {ionIndex, "b-" + labelBack};
                ERR_RETURN
            }
            else if (ionLabel.contains("^")) {
                ionLabelSplit = ionLabel.split("^", QString::SkipEmptyParts);
                e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

                const QString labelFront = ionLabelSplit.front().replace("b", "");

                const QString labelBack = ionLabelSplit.back();
                e = ErrorUtils::toInt(labelFront, &ionIndex); ree

                *ionInfo = {ionIndex, "b^" + labelBack};
                ERR_RETURN
            }
            else {

                QString ionLabelCopy = ionLabel;
                ionLabelCopy = ionLabelCopy.replace("b", "");
                e = ErrorUtils::toInt(ionLabelCopy, &ionIndex); ree
                *ionInfo = {ionIndex, "b"};
                ERR_RETURN
            }
        }

        else if (ionLabel.contains("a")) {
            QString ionLabelCopy = ionLabel;
            ionLabelCopy = ionLabelCopy.replace("a", "");
            e = ErrorUtils::toInt(ionLabelCopy, &ionIndex); ree
            *ionInfo = {ionIndex, "a"};
            ERR_RETURN
        }

        rrr(eValueError);
    }

    static void sortMS2IonsMzAsc(QVector<MS2Ion> *ms2Ions) {

        const auto sortMzAsc = [](const MS2Ion &l, const MS2Ion &r){
            return l.mz < r.mz;
        };

        std::sort(ms2Ions->begin(), ms2Ions->end(), sortMzAsc);
    }

    static void filterMS2IonsByMz(
            double mzMin,
            double mzMax,
            QVector<MS2Ion> *ms2Ions
            ) {

        const auto terminatorLogic = [mzMin, mzMax](const MS2Ion &ion){
            return !(mzMin < ion.mz && ion.mz <= mzMax);
        };

        const auto terminator = std::remove_if(
                ms2Ions->begin(),
                ms2Ions->end(),
                terminatorLogic
                );

        ms2Ions->erase(terminator, ms2Ions->end());
    }

    static ScanPoints ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions) {

        const auto convLog = [](const MS2Ion &ion){
            return QPointF(ion.mz, ion.intensity);
        };

        ScanPoints predPoints;
        std::transform(
                ms2Ions.begin(),
                ms2Ions.end(),
                std::back_inserter(predPoints),
                convLog
                );

        return predPoints;
    }

};

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

    const QString Y_IONS =  QStringLiteral("yIons");
    const QString B_IONS =  QStringLiteral("bIons");
    const QString Y2_IONS =  QStringLiteral("y2Ions");
    const QString B2_IONS =  QStringLiteral("b2Ions");
    const QString A_IONS =  QStringLiteral("aIons");
    const QString Y_NH3_IONS =  QStringLiteral("yNH3");
    const QString Y_H2O_IONS =  QStringLiteral("yH20");
    const QString B_NH3_IONS=  QStringLiteral("bNH3");
    const QString B_H2O_IONS =  QStringLiteral("bH20");
    const QString PRECURSOR_IONS =  QStringLiteral("precursorIons");

    static QString VERSION();
};


const extern UTILSLIB_EXPORTS GlobalSettings S_GLOBAL_SETTINGS;


#endif //GLOBALSETTINGS_H
