//
// Created by anichols on 4/3/23.
//

#ifndef PYTHIADIACPP_FRAGLIBREADER_H
#define PYTHIADIACPP_FRAGLIBREADER_H

#include <utility>

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

using namespace Error;

namespace FragLibReaderNamespace {

    const QString PEP_SEQ_CHRG_KEY = QStringLiteral("peptideSequenceChargeKey");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString MASS = QStringLiteral("mass");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString ION_LABLES = QStringLiteral("ionLabels");
    const QString IRT_LABEL = QStringLiteral("iRT");

    const QStringList keysToCheck = {
            PEP_SEQ_CHRG_KEY,
            MZ_VALS,
            INTENSITY_VALS,
            MASS,
            IS_DECOY,
            ION_LABLES,
            IRT_LABEL
    };
}

struct FILEREADERSLIB_EXPORTS FragLibReaderRow : public ParquetReaderInputBase {

    PeptideSequenceChargeKey peptideSequenceChargeKey;
    QVector<double> mzVals;
    QVector<double> intensityVals;
    QString ionLabels;
    double mass = -1.0;
    int isDecoy = 0;
    double iRT = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace FragLibReaderNamespace;

        return {
            {PEP_SEQ_CHRG_KEY, QVariant(peptideSequenceChargeKey)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {MASS, QVariant(mass)},
            {IS_DECOY, QVariant(isDecoy)},
            {ION_LABLES, QVariant(ionLabels)},
            {IRT_LABEL, QVariant(iRT)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace FragLibReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
                );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideSequenceChargeKey = dataMap.value(PEP_SEQ_CHRG_KEY).toString();
        mzVals = bytesArrayToQVector<double>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());
        mass = dataMap.value(MASS).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toInt();
        ionLabels = dataMap.value(ION_LABLES).toString();
        iRT = dataMap.value(IRT_LABEL).toDouble();

        ERR_RETURN
    }

};

struct CandidatePeptide {
    PeptideStringWithMods peptideStringWithMods;
    Charge charge = -1;
    QVector<MS2Ion> ms2Ions;
    bool isDecoy = false;
    double mass = -1.0;
    double iRt = -1.0;
    int totalFragmentCount;
};


class FILEREADERSLIB_EXPORTS FragLibReader {


public:

    FragLibReader();
    ~FragLibReader() = default;

    Err init(const QString &fragLibFilePath);

    int libarySize();

    Err getMS2Ions(QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide);

    Err getMS2IonsTopN(
            int topNMs2Ions,
            double mzMin,
            double mzMax,
            QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
    );

    Err getMS2IonsTopN(
            const QMap<Index, bool> &selectionList,
            int topNMs2Ions,
            double mzMin,
            double mzMax,
            QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
    );

    static Err getMS2Ions(
            const QString &fragLibFilePath,
            double massStart,
            double massEnd,
            double mzMin,
            double mzMax,
            int topNMs2Ions,
            QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
    );

    static void filterMs2IonsByMz(
            double mzStart,
            double mzEnd,
            QVector<MS2Ion> *ms2Ions
    );

    static void filterMs2IonsByIntensity(
            double intensityThreshold,
            QVector<MS2Ion> *ms2Ions
            );

    static void getTopNMostIntenseMs2Ions(
            int topNMs2Ions,
            double mzMin,
            double mzMax,
            QVector<MS2Ion> *ms2Ions
    );

    static Err buildFragIonLibForCandidates(
            const QString &fragLibUri,
            int chargeMin,
            int chargeMax,
            double mzTargetMin,
            double mzTargetMax,
            QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
            );

    static Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
    );

    static Err fragLibReaderRowsToMs2IonsMap(
            const QVector<FragLibReaderRow> &fragLibReaderRows,
            int topNMs2Ions,
            double mzMin,
            double mzMax,
            QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
    );

    static Err generateFragmentFrequencies(
            const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
            double ppmTol,
            QMap<MzHashed, FrequencyPercent> *fragmentFrequencies
    );

    static Err convertDIANNLibToFragLib(const QString &specLibFilePath);

    static Err mutateCandidatePeptideTarget(
            const CandidatePeptide &candidatePeptideTarget,
            CandidatePeptide *candidatePeptideDecoy
    );

private:

    QString m_fragLibFilePath;
    double m_mzMin;
    double m_mzMax;

    QVector<FragLibReaderRow> m_fragLibReaderRows;

};




#endif //PYTHIADIACPP_FRAGLIBREADER_H
