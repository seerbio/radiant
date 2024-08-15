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

    const QString SCAN_TIMES_1 = QStringLiteral("scanTimesMz1");
    const QString SCAN_TIMES_2 = QStringLiteral("scanTimesMz2");
    const QString SCAN_TIMES_3 = QStringLiteral("scanTimesMz3");
    const QString SCAN_TIMES_4 = QStringLiteral("scanTimesMz4");
    const QString SCAN_TIMES_5 = QStringLiteral("scanTimesMz5");
    const QString SCAN_TIMES_6 = QStringLiteral("scanTimesMz6");

    const QString INTENSITY_VALS_MZ1 = QStringLiteral("intensityValsMz1");
    const QString INTENSITY_VALS_MZ2 = QStringLiteral("intensityValsMz2");
    const QString INTENSITY_VALS_MZ3 = QStringLiteral("intensityValsMz3");
    const QString INTENSITY_VALS_MZ4 = QStringLiteral("intensityValsMz4");
    const QString INTENSITY_VALS_MZ5 = QStringLiteral("intensityValsMz5");
    const QString INTENSITY_VALS_MZ6 = QStringLiteral("intensityValsMz6");

    const QString CLASSIFIER_SCORE = QStringLiteral("classifierScore");
    const QString DISC_SCORE = QStringLiteral("discScore");
    const QString Q_VALUE = QStringLiteral("qValue");

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
        SCAN_TIMES_1,
        SCAN_TIMES_2,
        SCAN_TIMES_3,
        SCAN_TIMES_4,
        SCAN_TIMES_5,
        SCAN_TIMES_6,
        INTENSITY_VALS_MZ1,
        INTENSITY_VALS_MZ2,
        INTENSITY_VALS_MZ3,
        INTENSITY_VALS_MZ4,
        INTENSITY_VALS_MZ5,
        INTENSITY_VALS_MZ6,
        CLASSIFIER_SCORE,
        DISC_SCORE,
        Q_VALUE
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
    QVector<float> scanTimesMz1;
    QVector<float> intensityValsMz1;
    QVector<float> scanTimesMz2;
    QVector<float> intensityValsMz2;
    QVector<float> scanTimesMz3;
    QVector<float> intensityValsMz3;
    QVector<float> scanTimesMz4;
    QVector<float> intensityValsMz4;
    QVector<float> scanTimesMz5;
    QVector<float> intensityValsMz5;
    QVector<float> scanTimesMz6;
    QVector<float> intensityValsMz6;

    float classifierScore = -1.0;
    float discScore = -1.0;
    float qValue = -1.0;

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

            {SCAN_TIMES_1, QVariant(qVectorToQByteArray(scanTimesMz1))},
            {SCAN_TIMES_2, QVariant(qVectorToQByteArray(scanTimesMz2))},
            {SCAN_TIMES_3, QVariant(qVectorToQByteArray(scanTimesMz3))},
            {SCAN_TIMES_4, QVariant(qVectorToQByteArray(scanTimesMz4))},
            {SCAN_TIMES_5, QVariant(qVectorToQByteArray(scanTimesMz5))},
            {SCAN_TIMES_6, QVariant(qVectorToQByteArray(scanTimesMz6))},

            {INTENSITY_VALS_MZ1, QVariant(qVectorToQByteArray(intensityValsMz1))},
            {INTENSITY_VALS_MZ2, QVariant(qVectorToQByteArray(intensityValsMz2))},
            {INTENSITY_VALS_MZ3, QVariant(qVectorToQByteArray(intensityValsMz3))},
            {INTENSITY_VALS_MZ4, QVariant(qVectorToQByteArray(intensityValsMz4))},
            {INTENSITY_VALS_MZ5, QVariant(qVectorToQByteArray(intensityValsMz5))},
            {INTENSITY_VALS_MZ6, QVariant(qVectorToQByteArray(intensityValsMz6))},

            {CLASSIFIER_SCORE, QVariant(classifierScore)},
            {DISC_SCORE, QVariant(discScore)},
            {Q_VALUE, QVariant(qValue)},
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

        scanTimesMz1 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_1).toByteArray());
        scanTimesMz2 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_2).toByteArray());
        scanTimesMz3 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_3).toByteArray());
        scanTimesMz4 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_4).toByteArray());
        scanTimesMz5 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_5).toByteArray());
        scanTimesMz6 = bytesArrayToQVector<float>(dataMap.value(SCAN_TIMES_6).toByteArray());

        intensityValsMz1 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ1).toByteArray());
        intensityValsMz2 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ2).toByteArray());
        intensityValsMz3 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ3).toByteArray());
        intensityValsMz4 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ4).toByteArray());
        intensityValsMz5 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ5).toByteArray());
        intensityValsMz6 = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS_MZ6).toByteArray());

        classifierScore = dataMap.value(CLASSIFIER_SCORE).toFloat();
        discScore = dataMap.value(DISC_SCORE).toFloat();
        qValue = dataMap.value(Q_VALUE).toFloat();

        ERR_RETURN
    }

};


class FILEREADERSLIB_EXPORTS QuanReader {

};



#endif //QUANREADER_H
