//
// Created by anichols on 4/26/23.
//

#ifndef PYTHIADIACPP_PSMSREADER_H
#define PYTHIADIACPP_PSMSREADER_H

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

#include <QMap>
#include <vector>

using namespace Error;

namespace PSMsReaderRowNamespace {

    const QString FRAME_INDEX = QStringLiteral("frameIndex");
    const QString SCAN_NUMBER = QStringLiteral("scanNumber");
    const QString CHARGE = QStringLiteral("charge");
    const QString UNIQUE_TARGET_KEY = QStringLiteral("uniqueTargetKey");
    const QString PEPTIDE_W_MODS = QStringLiteral("peptideWithMods");
    const QString SCORE = QStringLiteral("score");
    const QString FRAME_RANK_SCORE = QStringLiteral("frameRankScore");
    const QString DISC_SCORE = QStringLiteral("discScore");
    const QString FRAME_RANK_DISC_SCORE = QStringLiteral("frameRankDiscScore");
    const QString IS_DECOY = QStringLiteral("isDecoy");

    const QStringList keysToCheck = {
            FRAME_INDEX,
            SCAN_NUMBER ,
            CHARGE,
            UNIQUE_TARGET_KEY,
            PEPTIDE_W_MODS,
            SCORE,
            FRAME_RANK_SCORE ,
            DISC_SCORE ,
            FRAME_RANK_DISC_SCORE,
            IS_DECOY
    };
}

struct FILEREADERSLIB_EXPORTS PSMsReaderRow : public ParquetReaderInputBase {

    FrameIndex frameIndex = -1;
    ScanNumber scanNumber = -1;
    Charge charge = -1;
    UniqueMsInfoScanKey uniqueMsInfoScanKey;

    PeptideStringWithMods peptideStringWithMods;

    Score score = -1.0;
    int frameRankScore = -1;

    DiscScore discScore = -1.0;
    int frameRankDiscScore = -1;

    bool isDecoy = false;

    QMap<QString, QVariant> map() override {

        using namespace PSMsReaderRowNamespace;

        return {
                {FRAME_INDEX, QVariant(frameIndex)},
                {SCAN_NUMBER, QVariant(scanNumber)},
                {CHARGE, QVariant(charge)},
                {UNIQUE_TARGET_KEY, QVariant(uniqueMsInfoScanKey)},
                {PEPTIDE_W_MODS, QVariant(peptideStringWithMods)},
                {SCORE, QVariant(score)},
                {FRAME_RANK_SCORE, QVariant(frameRankScore)},
                {DISC_SCORE, QVariant(discScore)},
                {FRAME_RANK_DISC_SCORE, QVariant(frameRankDiscScore)},
                {IS_DECOY, QVariant(isDecoy)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace PSMsReaderRowNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frameIndex = dataMap.value(FRAME_INDEX).toInt();
        scanNumber = dataMap.value(SCAN_NUMBER).toInt();
        charge = dataMap.value(CHARGE).toInt();
        uniqueMsInfoScanKey = dataMap.value(UNIQUE_TARGET_KEY).toString();
        peptideStringWithMods = dataMap.value(PEPTIDE_W_MODS).toString();
        score = dataMap.value(SCORE).toDouble();
        frameRankScore = dataMap.value(FRAME_RANK_SCORE).toInt();
        discScore = dataMap.value(DISC_SCORE).toDouble();
        frameRankDiscScore = dataMap.value(FRAME_RANK_DISC_SCORE).toInt();
        isDecoy = dataMap.value(IS_DECOY).toBool();

        ERR_RETURN
    }

};



class FILEREADERSLIB_EXPORTS PSMsReader {


public:

    PSMsReader() = default;
    ~PSMsReader() = default;





};


#endif //PYTHIADIACPP_PSMSREADER_H
