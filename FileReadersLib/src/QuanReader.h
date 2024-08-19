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

    const QString SCAN_TIMES = QStringLiteral("scanTimes");

    const QString INTENSITY_VALS_MZ1 = QStringLiteral("intensityValsMz1");
    const QString INTENSITY_VALS_MZ2 = QStringLiteral("intensityValsMz2");
    const QString INTENSITY_VALS_MZ3 = QStringLiteral("intensityValsMz3");
    const QString INTENSITY_VALS_MZ4 = QStringLiteral("intensityValsMz4");
    const QString INTENSITY_VALS_MZ5 = QStringLiteral("intensityValsMz5");
    const QString INTENSITY_VALS_MZ6 = QStringLiteral("intensityValsMz6");
    const QString INTENSITY_VALS_MZ7 = QStringLiteral("intensityValsMz7");
    const QString INTENSITY_VALS_MZ8 = QStringLiteral("intensityValsMz8");
    const QString INTENSITY_VALS_MZ9 = QStringLiteral("intensityValsMz9");
    const QString INTENSITY_VALS_MZ10 = QStringLiteral("intensityValsMz10");
    const QString INTENSITY_VALS_MZ11 = QStringLiteral("intensityValsMz11");
    const QString INTENSITY_VALS_MZ12 = QStringLiteral("intensityValsMz12");

    const QString CLASSIFIER_SCORE = QStringLiteral("classifierScore");
    const QString DISC_SCORE = QStringLiteral("discScore");
    const QString Q_VALUE = QStringLiteral("qValue");
    const QString MZ_INTERF = QStringLiteral("mzInterferences");

    const QString IS_DECOY = QStringLiteral("isDecoy");

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
        SCAN_TIMES,
        INTENSITY_VALS_MZ1,
        INTENSITY_VALS_MZ2,
        INTENSITY_VALS_MZ3,
        INTENSITY_VALS_MZ4,
        INTENSITY_VALS_MZ5,
        INTENSITY_VALS_MZ6,
        CLASSIFIER_SCORE,
        DISC_SCORE,
        Q_VALUE,
        MZ_INTERF,
        IS_DECOY
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
    QVector<float> scanTimes;
    QVector<float> intensityValsMz1;
    QVector<float> intensityValsMz2;
    QVector<float> intensityValsMz3;
    QVector<float> intensityValsMz4;
    QVector<float> intensityValsMz5;
    QVector<float> intensityValsMz6;
    QVector<float> intensityValsMz7;
    QVector<float> intensityValsMz8;
    QVector<float> intensityValsMz9;
    QVector<float> intensityValsMz10;
    QVector<float> intensityValsMz11;
    QVector<float> intensityValsMz12;

    QVector<float> mzInterferences;

    float classifierScore = -1.0;
    float discScore = -1.0;
    float qValue = -1.0;

    bool isDecoy = false;

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

            {SCAN_TIMES, QVariant(qVectorToQByteArray(scanTimes))},

            {INTENSITY_VALS_MZ1, QVariant(qVectorToQByteArray(intensityValsMz1))},
            {INTENSITY_VALS_MZ2, QVariant(qVectorToQByteArray(intensityValsMz2))},
            {INTENSITY_VALS_MZ3, QVariant(qVectorToQByteArray(intensityValsMz3))},
            {INTENSITY_VALS_MZ4, QVariant(qVectorToQByteArray(intensityValsMz4))},
            {INTENSITY_VALS_MZ5, QVariant(qVectorToQByteArray(intensityValsMz5))},
            {INTENSITY_VALS_MZ6, QVariant(qVectorToQByteArray(intensityValsMz6))},
            {INTENSITY_VALS_MZ7, QVariant(qVectorToQByteArray(intensityValsMz7))},
            {INTENSITY_VALS_MZ8, QVariant(qVectorToQByteArray(intensityValsMz8))},
            {INTENSITY_VALS_MZ9, QVariant(qVectorToQByteArray(intensityValsMz9))},
            {INTENSITY_VALS_MZ10, QVariant(qVectorToQByteArray(intensityValsMz10))},
            {INTENSITY_VALS_MZ11, QVariant(qVectorToQByteArray(intensityValsMz11))},
            {INTENSITY_VALS_MZ12, QVariant(qVectorToQByteArray(intensityValsMz12))},

            {MZ_INTERF, QVariant(qVectorToQByteArray(mzInterferences))},

            {CLASSIFIER_SCORE, QVariant(classifierScore)},
            {DISC_SCORE, QVariant(discScore)},
            {Q_VALUE, QVariant(qValue)},

            {IS_DECOY, QVariant(isDecoy)}
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

        scanTimes = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES).toByteArray());

        intensityValsMz1 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ1).toByteArray());
        intensityValsMz2 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ2).toByteArray());
        intensityValsMz3 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ3).toByteArray());
        intensityValsMz4 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ4).toByteArray());
        intensityValsMz5 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ5).toByteArray());
        intensityValsMz6 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ6).toByteArray());
        intensityValsMz7 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ7).toByteArray());
        intensityValsMz8 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ8).toByteArray());
        intensityValsMz9 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ9).toByteArray());
        intensityValsMz10 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ10).toByteArray());
        intensityValsMz11 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ11).toByteArray());
        intensityValsMz12 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ12).toByteArray());

        mzInterferences = bytesArrayToQVector<float>(dataMap.value(MZ_INTERF).toByteArray());

        classifierScore = dataMap.value(CLASSIFIER_SCORE).toFloat();
        discScore = dataMap.value(DISC_SCORE).toFloat();
        qValue = dataMap.value(Q_VALUE).toFloat();

        isDecoy = dataMap.value(IS_DECOY).toBool();

        ERR_RETURN
    }

};


class FILEREADERSLIB_EXPORTS QuanReader {

};



#endif //QUANREADER_H
