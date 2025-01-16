//
// Created by andrewnichols on 11/15/24.
//

#ifndef FRAGLIBREADERROW_H
#define FRAGLIBREADERROW_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "ParquetReader.h"

using namespace Error;

namespace FragLibReaderNamespace {

    const QString FRAG_LIB_FF_SUFFIX = QStringLiteral("fragLibFF");
    const QString FRAG_LIB_TSV_FF_SUFFIX = QStringLiteral("tsv");
    const QString FRAG_LIB_CSV_FF_SUFFIX = QStringLiteral("csv");
    const QString SPEC_LIB_CSV_FF_SUFFIX = QStringLiteral("speclib");

    const QString PEP_SEQ_CHRG_KEY = QStringLiteral("peptideSequenceChargeKey");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString MASS = QStringLiteral("mass");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString ION_LABLES = QStringLiteral("ionLabels");
    const QString IRT_LABEL = QStringLiteral("iRT");
    const QString IM_LABEL = QStringLiteral("iM");

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

/**
 * See ParquetReaderInputBase for documentation
 */
struct FILEREADERSLIB_EXPORTS FragLibReaderRow : public ParquetReaderInputBase {

    PeptideSequenceChargeKey peptideSequenceChargeKey;
    QVector<float> mzVals;
    QVector<float> intensityVals;
    QString ionLabels;
    double mass = -1.0;
    int isDecoy = 0;
    double iRT = -1.0;
    double iM = -1.0;
    int precursorCharge = -1; //TODO include charge as column in fragLib files.
    QString proteinGroups;

    QMap<QString, QVariant> map() override {

        using namespace FragLibReaderNamespace;

        return {
            {PEP_SEQ_CHRG_KEY, QVariant(peptideSequenceChargeKey)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {MASS, QVariant(mass)},
            {IS_DECOY, QVariant(isDecoy)},
            {ION_LABLES, QVariant(ionLabels)},
            {IRT_LABEL, QVariant(iRT)},
            {IM_LABEL, QVariant(iM)}
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
        mzVals = bytesArrayToQVector<float>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS).toByteArray());
        mass = dataMap.value(MASS).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toInt();
        ionLabels = dataMap.value(ION_LABLES).toString();
        iRT = dataMap.value(IRT_LABEL).toDouble();
        iM = dataMap.value(IM_LABEL).toDouble();

        const QStringList peptideSequenceChargeKeySplit
            = peptideSequenceChargeKey.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP);

        //TODO after including charge as column in fragLib, delete this code.
        constexpr int expectedSplitSize = 2;
        e = ErrorUtils::isEqual(peptideSequenceChargeKeySplit.size(), expectedSplitSize); ree;
        e = ErrorUtils::toInt(peptideSequenceChargeKeySplit.at(1), &precursorCharge); ree;

        ERR_RETURN
    }

};





#endif //FRAGLIBREADERROW_H
