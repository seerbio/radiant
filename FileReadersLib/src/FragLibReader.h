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
    QVector<float> intensityVals;
    QString ionLabels;
    double mass = -1.0;
    int isDecoy = 0;
    double iRT = -1.0;
    int precursorCharge = -1; //TODO include charge as column in fragLib files.

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
        intensityVals = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS).toByteArray());
        mass = dataMap.value(MASS).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toInt();
        ionLabels = dataMap.value(ION_LABLES).toString();
        iRT = dataMap.value(IRT_LABEL).toDouble();

        const QStringList peptideSequenceChargeKeySplit
            = peptideSequenceChargeKey.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP);

        //TODO after including charge as column in fragLib, delete this code.
        const int expectedSplitSize = 2;
        e = ErrorUtils::isEqual(peptideSequenceChargeKeySplit.size(), expectedSplitSize); ree;
        e = ErrorUtils::toInt(peptideSequenceChargeKeySplit.at(1), &precursorCharge); ree;

        ERR_RETURN
    }

};


class FILEREADERSLIB_EXPORTS FragLibReader {

public:

    FragLibReader() = default;
    ~FragLibReader() = default;

    static Err getFragLibReaderRows(
            const QString &fragLibFilePath,
            double massStart,
            double massEnd,
            QVector<FragLibReaderRow> *fragLibReaderRows
    );

};




#endif //PYTHIADIACPP_FRAGLIBREADER_H
