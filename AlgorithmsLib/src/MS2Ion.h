//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_MS2ION_H
#define PYTHIADIACPP_MS2ION_H


#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class ALGORITHMSLIB_EXPORTS MS2Ion {

public:

    double mz = -1.0;
    double intensity = -1.0;
    IRT iRT = -1.0;
    QString ionLabel;
    int rank = -1;
    int charge = -1;

    MS2Ion() = default;
    ~MS2Ion() = default;

    friend QDebug operator<<(QDebug dbg, const MS2Ion& obj) {
        dbg.nospace() << "MS2Ion(" << obj.mz << ", " << obj.intensity << ") ";
        return dbg;
    }

    //TODO write unit tests
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

    static void sortMS2IonsIntensityDesc(QVector<MS2Ion> *ms2Ions) {

        const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
            return l.intensity < r.intensity;
        };

        std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);
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

        const auto convLogic = [](const MS2Ion &ion){
            return QPointF(ion.mz, ion.intensity);
        };

        ScanPoints predPoints;
        std::transform(
                ms2Ions.begin(),
                ms2Ions.end(),
                std::back_inserter(predPoints),
                convLogic
        );

        return predPoints;
    }

};


#endif //PYTHIADIACPP_MS2ION_H
