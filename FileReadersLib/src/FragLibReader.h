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

    const QStringList keysToCheck = {
            PEP_SEQ_CHRG_KEY,
            MZ_VALS,
            INTENSITY_VALS,
            MASS,
            IS_DECOY,
            ION_LABLES
    };
}

struct FILEREADERSLIB_EXPORTS FragLibReaderRow : public ParquetReaderInputBase {

    PeptideSequenceChargeKey peptideSequenceChargeKey;
    QVector<double> mzVals;
    QVector<double> intensityVals;
    QString ionLabels;
    double mass = -1.0;
    int isDecoy = 0;

    QMap<QString, QVariant> map() override {

        using namespace FragLibReaderNamespace;

        return {
            {PEP_SEQ_CHRG_KEY, QVariant(peptideSequenceChargeKey)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {MASS, QVariant(mass)},
            {IS_DECOY, QVariant(isDecoy)},
            {ION_LABLES, QVariant(ionLabels)}
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

        ERR_RETURN
    }

};

struct MS2IonsSeparated {

public:

    QMap<IonIndex, MS2Ion> yIons;
    QMap<IonIndex, MS2Ion> bIons;
    QMap<IonIndex, MS2Ion> y2Ions;
    QMap<IonIndex, MS2Ion> b2Ions;
    QMap<IonIndex, MS2Ion> aIons;
    QMap<IonIndex, MS2Ion> yNH3Ions;
    QMap<IonIndex, MS2Ion> yH2OIons;
    QMap<IonIndex, MS2Ion> bNH3Ions;
    QMap<IonIndex, MS2Ion> bH2OIons;
    QMap<IonIndex, MS2Ion> precursorIons;

    [[nodiscard]] QVector<MS2Ion> getTopXMS2Ions(
            int topX,
            int ionIndexCutoffThreshold
            ) const {

        QVector<MS2Ion> allMS2Ions;

        const QVector<MS2Ion> yIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                yIons,
                ionIndexCutoffThreshold
                );
        allMS2Ions.append(yIonsThresholded);

        const QVector<MS2Ion> y2IonsThresholded = returnMS2IonAboveIonIndexCutoff(
                y2Ions,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(y2IonsThresholded);

        const QVector<MS2Ion> yNH3IonsThresholded = returnMS2IonAboveIonIndexCutoff(
                yNH3Ions,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(yNH3IonsThresholded);

        const QVector<MS2Ion> yH2OIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                yH2OIons,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(yH2OIonsThresholded);

        const QVector<MS2Ion> bIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                bIons,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(bIonsThresholded);

        const QVector<MS2Ion> b2IonsThresholded = returnMS2IonAboveIonIndexCutoff(
                b2Ions,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(b2IonsThresholded);

        const QVector<MS2Ion> bNH3IonsThresholded = returnMS2IonAboveIonIndexCutoff(
                bNH3Ions,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(bNH3IonsThresholded);

        const QVector<MS2Ion> bH2OIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                bH2OIons,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(bH2OIonsThresholded);

        const QVector<MS2Ion> aIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                aIons,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(aIonsThresholded);

        const QVector<MS2Ion> precursorIonsThresholded = returnMS2IonAboveIonIndexCutoff(
                precursorIons,
                ionIndexCutoffThreshold
        );
        allMS2Ions.append(precursorIonsThresholded);

        std::sort(
                allMS2Ions.rbegin(),
                allMS2Ions.rend(),
                [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;}
                );

        topX = std::min(topX, allMS2Ions.size());

        allMS2Ions.resize(topX);
        return allMS2Ions;
    }

private:

    static QVector<MS2Ion> returnMS2IonAboveIonIndexCutoff(
            const QMap<IonIndex, MS2Ion> &ions,
            int ionIndexCutoffThreshold
            ) {

        QVector<MS2Ion> ms2IonsThresholded;

        for (auto it = ions.begin(); it != ions.end(); it++) {

            const IonIndex ionIndex = it.key();

            if (ionIndex < ionIndexCutoffThreshold) {
                continue;
            }

            ms2IonsThresholded.push_back(it.value());
        }

        return ms2IonsThresholded;
    }


};


class FILEREADERSLIB_EXPORTS FragLibReader {


public:

    FragLibReader() = default;
    ~FragLibReader() = default;

    Err init(const QString &fragLibFilePath);

    Err getMS2Ions(
            double massStart,
            double massEnd,
            QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideSequenceChargeKeyVsMS2Ions,
            QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy
    );

    Err getMS2Ions(
            double massStart,
            double massEnd,
            QMap<PeptideSequenceChargeKey, MS2IonsSeparated> *peptideSequenceChargeKeyVsMS2Ions,
            QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsMS2IonsSeparated
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
            QVector<MS2Ion> *ms2Ions
    );

    static Err buildFragIonLibForCandidates(
            const QString &fragLibUri,
            int chargeMin,
            int chargeMax,
            double mzTargetMin,
            double mzTargetMax,
            QMap<PeptideStringWithMods, MS2IonsSeparated> *fragPreds,
            QMap<PeptideStringWithMods, bool> *fragPredsIsDecoy
            );


private:

    QString m_fragLibFilePath;

};




#endif //PYTHIADIACPP_FRAGLIBREADER_H
