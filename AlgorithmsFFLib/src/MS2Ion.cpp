//
// Created by anichols on 10/17/23.
//

#include "MS2Ion.h"

Err MS2Ion::getIonLabelInfo(QPair<IonIndex, IonType> *ionInfo) const {

    ERR_INIT

    const int expectedSplitSize = 2;

    IonIndex ionIndex;

    if (ionLabel.contains("p")) {
        *ionInfo = {0, ionLabel};
        ERR_RETURN
    }
    else if (ionLabel.contains("y")) {

        QStringList ionLabelSplit;
        if (ionLabel.contains("-")) {

            ionLabelSplit = ionLabel.split("-", Qt::SkipEmptyParts);
            e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

            const QString labelFront = ionLabelSplit.front().replace("y", "");

            const QString labelBack = ionLabelSplit.back();
            e = ErrorUtils::toInt(labelFront, &ionIndex); ree

            *ionInfo = {ionIndex, "y-" + labelBack};
            ERR_RETURN
        }
        else if (ionLabel.contains("^")) {
            ionLabelSplit = ionLabel.split("^", Qt::SkipEmptyParts);
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

            ionLabelSplit = ionLabel.split("-", Qt::SkipEmptyParts);
            e = ErrorUtils::isEqual(ionLabelSplit.size(), expectedSplitSize); ree;

            const QString labelFront = ionLabelSplit.front().replace("b", "");

            const QString labelBack = ionLabelSplit.back();
            e = ErrorUtils::toInt(labelFront, &ionIndex); ree

            *ionInfo = {ionIndex, "b-" + labelBack};
            ERR_RETURN
        }
        else if (ionLabel.contains("^")) {
            ionLabelSplit = ionLabel.split("^", Qt::SkipEmptyParts);
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

void MS2Ion::sortMS2IonsMzAsc(QVector<MS2Ion> *ms2Ions) {

    const auto sortMzAsc = [](const MS2Ion &l, const MS2Ion &r){
        return l.mz < r.mz;
    };

    std::sort(ms2Ions->begin(), ms2Ions->end(), sortMzAsc);
}

void MS2Ion::sortMS2IonsIntensityDesc(QVector<MS2Ion> *ms2Ions) {

    const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
        return l.intensity < r.intensity;
    };

    std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);
}

void MS2Ion::filterMS2IonsByMz(
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

ScanPoints MS2Ion::ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions) {

    const auto convLogic = [](const MS2Ion &ion){
        return ScanPoint (ion.mz, ion.intensity);
    };

    ScanPoints predPoints;
    predPoints.reserve(ms2Ions.size());
    std::transform(
            ms2Ions.begin(),
            ms2Ions.end(),
            std::back_inserter(predPoints),
            convLogic
    );

    return predPoints;
}
