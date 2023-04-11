//
// Created by anichols on 4/10/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
#define PYTHIADIACPP_FRAGLIBRARYTRONDIA_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "PeptideMassRTree.h"
#include "PythiaParameterReader.h"


struct ALGORITHMSLIB_EXPORTS MS2Ion {
    double mz = -1.0;
    double intensity = -1.0;
    QString ionLabel;

    friend QDebug &operator <<(QDebug &d, const MS2Ion &ms2Ion) {
        d << "mz: " << ms2Ion.mz << ", ";
        d << "intz: " << ms2Ion.intensity << ", ";
        d << "ion: " << ms2Ion.ionLabel;
        return d;
    }
};


class ALGORITHMSLIB_EXPORTS FragLibraryTronDIA {

public:

    friend class FragLibraryTronDIATests;

    FragLibraryTronDIA() = default;
    ~FragLibraryTronDIA() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibFilePath
            );

    Err getMS2Ions(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            QVector<MS2Ion> *ms2Ions
            );

    Err getMS2Ions(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            QPair<double, double> mzStartMzEnd,
            QVector<MS2Ion> *ms2Ions
            );

    Err getMS2Ions(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            QPair<double, double> mzStartMzEnd,
            int topNIntense,
            QVector<MS2Ion> *ms2Ions
            );

    Err getMS2Ions(
            double mzTargetStart,
            double mzTargetEnd,
            int topNIntense,
            QHash<PeptideStringWithMods, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions
            );

private:

    Err readFragLibFile(const QString &fragLibFilePath);

    void buildChargeVsIonLabels();


private:

    QHash<PeptideSequenceChargeKey, QVector<MS2Ion>> m_pepSeqChrgKeyVsMS2Ions;

    PythiaParameters m_params;

    QMap<Charge, IonLabels> m_chargeVsIonLabels;

    PeptideMassRTree m_peptideMassRTree;

};


#endif //PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
