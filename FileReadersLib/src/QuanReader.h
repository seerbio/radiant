//
// Created by andrewnichols on 8/13/24.
//

#ifndef QUANREADER_H
#define QUANREADER_H

#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

#include <QMap>
#include <vector>

using namespace Error;

namespace QuanReaderNamespace {

    const QString PEP_STR_W_MODS = QStringLiteral("PeptideStringWithMods");
    const QString CHARGE = QStringLiteral("Charge");
    const QString TARG_KEY = QStringLiteral("TargetKey");
    const QString MZ_SEARCHED_1 = QStringLiteral("MzSearched1");
    const QString MZ_SEARCHED_2 = QStringLiteral("MzSearched2");
    const QString MZ_SEARCHED_3 = QStringLiteral("MzSearched3");
    const QString MZ_SEARCHED_4 = QStringLiteral("MzSearched4");
    const QString MZ_SEARCHED_5 = QStringLiteral("MzSearched5");
    const QString MZ_SEARCHED_6 = QStringLiteral("MzSearched6");
    const QString MZ_SEARCHED_7 = QStringLiteral("MzSearched7");
    const QString MZ_SEARCHED_8 = QStringLiteral("MzSearched8");
    const QString MZ_SEARCHED_9 = QStringLiteral("MzSearched9");
    const QString MZ_SEARCHED_10 = QStringLiteral("MzSearched10");
    const QString MZ_SEARCHED_11 = QStringLiteral("MzSearched11");
    const QString MZ_SEARCHED_12 = QStringLiteral("MzSearched12");

    const QString SCAN_TIME = QStringLiteral("scanTime");
    const QString SCAN_TIME_START = QStringLiteral("scanTimeStart");
    const QString SCAN_TIME_END = QStringLiteral("scanTimeEnd");

    const QString INTENSITY_VAL_MZ1 = QStringLiteral("intensityValMz1");
    const QString INTENSITY_VAL_MZ2 = QStringLiteral("intensityValMz2");
    const QString INTENSITY_VAL_MZ3 = QStringLiteral("intensityValMz3");
    const QString INTENSITY_VAL_MZ4 = QStringLiteral("intensityValMz4");
    const QString INTENSITY_VAL_MZ5 = QStringLiteral("intensityValMz5");
    const QString INTENSITY_VAL_MZ6 = QStringLiteral("intensityValMz6");

    const QString COS_SIM_SUM_100 = QStringLiteral("cosineSimSum100");

    const QString CLASSIFIER_SCORE = QStringLiteral("classifierScore");
    const QString DISC_SCORE = QStringLiteral("discScore");
    const QString Q_VALUE = QStringLiteral("qValue");
    const QString MZ_INTERF = QStringLiteral("mzInterferences");

    const QString IS_DECOY = QStringLiteral("isDecoy");

    const QString RAW_INTENSITY = QStringLiteral("rawIntensity");

    const QString COS_SIM_SUM_100_ALT = QStringLiteral("cosineSimSum100Alt");
    const QString SCAN_TIME_PRED = QStringLiteral("scanTimePredicted");

    const QStringList keysToCheck = {
        PEP_STR_W_MODS,
        CHARGE,
        TARG_KEY,
        MZ_SEARCHED_1,
        MZ_SEARCHED_2,
        MZ_SEARCHED_3,
        MZ_SEARCHED_4,
        MZ_SEARCHED_5,
        MZ_SEARCHED_6,
        MZ_SEARCHED_7,
        MZ_SEARCHED_8,
        MZ_SEARCHED_9,
        MZ_SEARCHED_10,
        MZ_SEARCHED_11,
        MZ_SEARCHED_12,
        SCAN_TIME,
        SCAN_TIME_START,
        SCAN_TIME_END,
        INTENSITY_VAL_MZ1,
        INTENSITY_VAL_MZ2,
        INTENSITY_VAL_MZ3,
        INTENSITY_VAL_MZ4,
        INTENSITY_VAL_MZ5,
        INTENSITY_VAL_MZ6,
        RAW_INTENSITY,
        COS_SIM_SUM_100,
        CLASSIFIER_SCORE,
        DISC_SCORE,
        Q_VALUE,
        MZ_INTERF,
        IS_DECOY,
        COS_SIM_SUM_100_ALT,
        SCAN_TIME_PRED
    };

}//namespace

struct FILEREADERSLIB_EXPORTS QuanReaderRow : public ParquetReaderInputBase {

    QString peptideStringWithMods;
    int charge = -1;
    QString targetKey;
    float mzSearched1 = -1.0;
    float mzSearched2 = -1.0;
    float mzSearched3 = -1.0;
    float mzSearched4 = -1.0;
    float mzSearched5 = -1.0;
    float mzSearched6 = -1.0;
    float mzSearched7 = -1.0;
    float mzSearched8 = -1.0;
    float mzSearched9 = -1.0;
    float mzSearched10 = -1.0;
    float mzSearched11 = -1.0;
    float mzSearched12 = -1.0;

    float scanTime = -1.0;
    float scanTimeStart = -1.0;
    float scanTimeEnd = -1.0;
    float intensityValMz1 = -1.0;
    float intensityValMz2 = -1.0;
    float intensityValMz3 = -1.0;
    float intensityValMz4 = -1.0;
    float intensityValMz5 = -1.0;
    float intensityValMz6 = -1.0;

    float cosineSimSum100 = -1.0;

    float classifierScore = -1.0;
    float discScore = -1.0;
    float qValue = -1.0;

    float rawIntensity = -1.0;

    bool isDecoy = false;

    float scanTimePredicted = -1.0;


    QMap<QString, QVariant> map() override {

        using namespace QuanReaderNamespace;

        return {
            {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
            {CHARGE, QVariant(charge)},
            {TARG_KEY, QVariant(targetKey)},
            {MZ_SEARCHED_1, QVariant(mzSearched1)},
            {MZ_SEARCHED_2, QVariant(mzSearched2)},
            {MZ_SEARCHED_3, QVariant(mzSearched3)},
            {MZ_SEARCHED_4, QVariant(mzSearched4)},
            {MZ_SEARCHED_5, QVariant(mzSearched5)},
            {MZ_SEARCHED_6, QVariant(mzSearched6)},
            {MZ_SEARCHED_7, QVariant(mzSearched7)},
            {MZ_SEARCHED_8, QVariant(mzSearched8)},
            {MZ_SEARCHED_9, QVariant(mzSearched9)},
            {MZ_SEARCHED_10, QVariant(mzSearched10)},
            {MZ_SEARCHED_11, QVariant(mzSearched11)},
            {MZ_SEARCHED_12, QVariant(mzSearched12)},

            {SCAN_TIME, QVariant(scanTime)},
            {SCAN_TIME_START, QVariant(scanTimeStart)},
            {SCAN_TIME_END, QVariant(scanTimeEnd)},
            {INTENSITY_VAL_MZ1, QVariant(intensityValMz1)},
            {INTENSITY_VAL_MZ2, QVariant(intensityValMz2)},
            {INTENSITY_VAL_MZ3, QVariant(intensityValMz3)},
            {INTENSITY_VAL_MZ4, QVariant(intensityValMz4)},
            {INTENSITY_VAL_MZ5, QVariant(intensityValMz5)},
            {INTENSITY_VAL_MZ6, QVariant(intensityValMz6)},

            {COS_SIM_SUM_100, QVariant(cosineSimSum100)},

            {RAW_INTENSITY, QVariant(rawIntensity)},

            {CLASSIFIER_SCORE, QVariant(classifierScore)},
            {DISC_SCORE, QVariant(discScore)},
            {Q_VALUE, QVariant(qValue)},

            {IS_DECOY, QVariant(isDecoy)},

            {SCAN_TIME_PRED, QVariant(scanTimePredicted)},

        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace QuanReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        targetKey = dataMap.value(TARG_KEY).toString();
        mzSearched1 = dataMap.value(MZ_SEARCHED_1).toFloat();
        mzSearched2 = dataMap.value(MZ_SEARCHED_2).toFloat();
        mzSearched3 = dataMap.value(MZ_SEARCHED_3).toFloat();
        mzSearched4 = dataMap.value(MZ_SEARCHED_4).toFloat();
        mzSearched5 = dataMap.value(MZ_SEARCHED_5).toFloat();
        mzSearched6 = dataMap.value(MZ_SEARCHED_6).toFloat();
        mzSearched7 = dataMap.value(MZ_SEARCHED_7).toFloat();
        mzSearched8 = dataMap.value(MZ_SEARCHED_8).toFloat();
        mzSearched9 = dataMap.value(MZ_SEARCHED_9).toFloat();
        mzSearched10 = dataMap.value(MZ_SEARCHED_10).toFloat();
        mzSearched11 = dataMap.value(MZ_SEARCHED_11).toFloat();
        mzSearched12 = dataMap.value(MZ_SEARCHED_12).toFloat();

        scanTime = dataMap.value(SCAN_TIME).toFloat();
        scanTimeStart = dataMap.value(SCAN_TIME_START).toFloat();
        scanTimeEnd = dataMap.value(SCAN_TIME_END).toFloat();
        intensityValMz1 = dataMap.value(INTENSITY_VAL_MZ1).toFloat();
        intensityValMz2 = dataMap.value(INTENSITY_VAL_MZ2).toFloat();
        intensityValMz3 = dataMap.value(INTENSITY_VAL_MZ3).toFloat();
        intensityValMz4 = dataMap.value(INTENSITY_VAL_MZ4).toFloat();
        intensityValMz5 = dataMap.value(INTENSITY_VAL_MZ5).toFloat();
        intensityValMz6 = dataMap.value(INTENSITY_VAL_MZ6).toFloat();

        cosineSimSum100 = dataMap.value(COS_SIM_SUM_100).toFloat();

        rawIntensity = dataMap.value(RAW_INTENSITY).toFloat();

        classifierScore = dataMap.value(CLASSIFIER_SCORE).toFloat();
        discScore = dataMap.value(DISC_SCORE).toFloat();
        qValue = dataMap.value(Q_VALUE).toFloat();

        isDecoy = dataMap.value(IS_DECOY).toBool();

        scanTimePredicted = dataMap.value(SCAN_TIME_PRED).toFloat();

        ERR_RETURN
    }

};


class FILEREADERSLIB_EXPORTS QuanReader {

};



#endif //QUANREADER_H
