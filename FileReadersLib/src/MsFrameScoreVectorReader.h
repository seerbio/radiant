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

    const QStringList keysToCheck = {
            COSINE_SIM_VEC,
            INTZ_FRM_IND,
            FOUNDS_FRM_IND,
            SCORE_FRM_IND,
            PEP_STR_W_MODS
    };
}

struct FILEREADERSLIB_EXPORTS MsFrameScoreVectorReaderRow : public ParquetReaderInputBase {

    QVector<double> cosineSimPerFrameIndexOfTargetVec;
    QVector<double> intensityPerFrameIndexOfTargetVec;
    QVector<int> foundIonsPerFrameIndexOfTargetVec;
    QVector<double> scorePerFrameIndexOfTargetVec;
    PeptideStringWithMods peptideStringWithMods;

    QMap<QString, QVariant> map() override {

        using namespace MsFrameScoreVectorReaderNamespace;

        return {
            {COSINE_SIM_VEC, QVariant(qVectorToQByteArray(cosineSimPerFrameIndexOfTargetVec))},
            {INTZ_FRM_IND, QVariant(qVectorToQByteArray(intensityPerFrameIndexOfTargetVec))},
            {FOUNDS_FRM_IND, QVariant(qVectorToQByteArray(foundIonsPerFrameIndexOfTargetVec))},
            {SCORE_FRM_IND, QVariant(qVectorToQByteArray(scorePerFrameIndexOfTargetVec))},
            {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
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

        ERR_RETURN
    }
};



#endif //PYTHIADIACPP_MSFRAMESCOREVECTORREADER_H
