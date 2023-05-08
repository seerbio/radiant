//
// Created by anichols on 4/12/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCOREVECTORREADER_H
#define PYTHIADIACPP_MSFRAMESCOREVECTORREADER_H

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"


using namespace Error;


namespace MsFrameScoreVectorReaderNamespace {

    const QString COSINE_SIM_VEC = QStringLiteral("cosineSimPerFrameIndexOfTargetVec");
    const QString INTZ_FRM_IND = QStringLiteral("intensityPerFrameIndexOfTargetVec");
    const QString FOUNDS_FRM_IND = QStringLiteral("foundIonsPerFrameIndexOfTargetVec");
    const QString SCORE_FRM_IND = QStringLiteral("scorePerFrameIndexOfTargetVec");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString SCORE_PEAK_START  = QStringLiteral("scorePeakStart");
    const QString SCORE_PEAK_END  = QStringLiteral("scorePeakEnd");
    const QString FRAME_INDEX_MAX_SCORE = QStringLiteral("frameIndexMaxScore");
    const QString CHARGE = QStringLiteral("charge");

    const QStringList keysToCheck = {
            COSINE_SIM_VEC,
            INTZ_FRM_IND,
            FOUNDS_FRM_IND,
            SCORE_FRM_IND,
            PEP_STR_W_MODS,
            SCORE_PEAK_START,
            SCORE_PEAK_END,
            FRAME_INDEX_MAX_SCORE,
            CHARGE
    };
}

struct FILEREADERSLIB_EXPORTS MsFrameScoreVectorReaderRow : public ParquetReaderInputBase {

    QVector<double> cosineSimPerFrameIndexOfTargetVec;
    QVector<double> intensityPerFrameIndexOfTargetVec;
    QVector<int> foundIonsPerFrameIndexOfTargetVec;
    QVector<double> scorePerFrameIndexOfTargetVec;
    PeptideStringWithMods peptideStringWithMods;
    int scorePeakStart = -1;
    int scorePeakEnd = -1;
    int frameIndexMaxScore = -1;
    int charge = -1;

    QMap<QString, QVariant> map() override {

        using namespace MsFrameScoreVectorReaderNamespace;

        return {
            {COSINE_SIM_VEC, QVariant(qVectorToQByteArray(cosineSimPerFrameIndexOfTargetVec))},
            {INTZ_FRM_IND, QVariant(qVectorToQByteArray(intensityPerFrameIndexOfTargetVec))},
            {FOUNDS_FRM_IND, QVariant(qVectorToQByteArray(foundIonsPerFrameIndexOfTargetVec))},
            {SCORE_FRM_IND, QVariant(qVectorToQByteArray(scorePerFrameIndexOfTargetVec))},
            {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
            {SCORE_PEAK_START, QVariant(scorePeakStart)},
            {SCORE_PEAK_END, QVariant(scorePeakEnd)},
            {FRAME_INDEX_MAX_SCORE, QVariant(frameIndexMaxScore)},
            {CHARGE, QVariant(charge)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MsFrameScoreVectorReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        cosineSimPerFrameIndexOfTargetVec = bytesArrayToQVector<double>(dataMap.value(COSINE_SIM_VEC).toByteArray());
        intensityPerFrameIndexOfTargetVec = bytesArrayToQVector<double>(dataMap.value(COSINE_SIM_VEC).toByteArray());
        foundIonsPerFrameIndexOfTargetVec = bytesArrayToQVector<int>(dataMap.value(FOUNDS_FRM_IND).toByteArray());
        scorePerFrameIndexOfTargetVec = bytesArrayToQVector<double>(dataMap.value(SCORE_FRM_IND).toByteArray());
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        scorePeakStart = dataMap.value(SCORE_PEAK_START).toInt();
        scorePeakEnd = dataMap.value(SCORE_PEAK_END).toInt();
        frameIndexMaxScore = dataMap.value(FRAME_INDEX_MAX_SCORE).toInt();
        charge = dataMap.value(CHARGE).toInt();

        ERR_RETURN
    }
};


#endif //PYTHIADIACPP_MSFRAMESCOREVECTORREADER_H
