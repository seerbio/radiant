//
// Created by anichols on 4/10/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
#define PYTHIADIACPP_FRAGLIBRARYTRONDIA_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "PythiaParameterReader.h"

class FragLibReaderRow;


struct ALGORITHMSLIB_EXPORTS MS2Ion {

    double mz = -1.0;
    double intensity = -1.0;

    MS2Ion() = default;

    MS2Ion(double mz, double intensity)
    : mz(mz) , intensity(intensity){}

    friend QDebug &operator <<(QDebug &d, const MS2Ion &ms2Ion) {
        d << "MS2Ion(" << ms2Ion.mz << "," << ms2Ion.intensity << ")";
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
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions
            );

    Err peptideStringWithModsIsDecoy(
            const PeptideStringWithMods &peptideStringWithMods,
            bool *isDecoy
            );

    bool isInit();

    [[nodiscard]] static ScanPoints ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions);

    Err buildTargetCandidatesForFrame(
            const QVector<QPair<double, double>> &percursorMzWindows,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *framePredictions
    );

private:

    Err readFragLibFile(const QString &fragLibFilePath);

    Err buildPeptideWithModsVsisDecoy(const QVector<FragLibReaderRow> &tandemLibraryRows);

private:

    QHash<PeptideSequenceChargeKey, QVector<MS2Ion>> m_pepSeqChrgKeyVsMS2Ions;
    QHash<PeptideStringWithMods, bool> m_peptideWithModsVsisDecoy;

    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
