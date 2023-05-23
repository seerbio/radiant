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
    const QString RESCORE = QStringLiteral("rescore");
    const QString T_TEST = QStringLiteral("tTest");
    const QString P_VAL = QStringLiteral("pVal");
    const QString FRAME_F_STAT = QStringLiteral("frameFStat");
    const QString P_VAL_FRAME_F_TEST = QStringLiteral("pValFrameFtest");
    const QString FRAME_ERROR = QStringLiteral("frameError");
    const QString MISSED_CLEAVAGES = QStringLiteral("missedCleavages");

    const QString SCORE_MEAN = QStringLiteral("scoreMean");
    const QString SCORE_MEDIAN = QStringLiteral("scoreMedian");
    const QString SCORE_STDEV = QStringLiteral("scoreStDev");
    const QString SCORE_MIN = QStringLiteral("scoreMin");
    const QString SCORE_MAX = QStringLiteral("scoreMax");
    const QString DISC_SCORE_MEAN = QStringLiteral("discScoreMean");
    const QString DISC_SCORE_MEDIAN = QStringLiteral("discScoreMedian");
    const QString DISC_SCORE_STDEV = QStringLiteral("discScoreStDev");
    const QString DISC_SCORE_MIN = QStringLiteral("discScoreMin");
    const QString DISC_SCORE_MAX = QStringLiteral("discScoreMax");
    const QString COSINE_SIM = QStringLiteral("cosineSim");
    const QString KL_DIV = QStringLiteral("klDiv");
    const QString FRACTION_FOUND = QStringLiteral("fractionFound");

    const QString IONS_FOUND = QStringLiteral("ionsFound");
    const QString FRAME_CAND_COUNT = QStringLiteral("frameCandidateCount");
    const QString PEPTIDE_SIZE = QStringLiteral("peptideSize");
    const QString MZ = QStringLiteral("mz");
    const QString MASS_THEO = QStringLiteral("massTheoretical");
    const QString MASS_FOUND = QStringLiteral("massFound");
    const QString MS1COSINE_SIM = QStringLiteral("ms1CosineSim");
    const QString PPM_DIFF_MS1 = QStringLiteral("ppmDiffMs1");
    const QString MS1_INTENSITY = QStringLiteral("ms1Intensity");
    const QString MONO_ISO_OFFSET = QStringLiteral("monoIsoOffset");
    const QString ISOTOPE_FOUND_CNT = QStringLiteral("isotopeFoundCount");
    const QString MZ_FOUND = QStringLiteral("mzFound");


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
            IS_DECOY,
            RESCORE,
            P_VAL,
            T_TEST,
            FRAME_F_STAT,
            P_VAL_FRAME_F_TEST,
            FRAME_ERROR,
            SCORE_MEAN,
            SCORE_MEDIAN,
            SCORE_STDEV,
            SCORE_MIN,
            SCORE_MAX,
            DISC_SCORE_MEAN,
            DISC_SCORE_MEDIAN,
            DISC_SCORE_STDEV,
            DISC_SCORE_MIN,
            DISC_SCORE_MAX,
            COSINE_SIM,
            KL_DIV,
            FRACTION_FOUND,
            IONS_FOUND,
            FRAME_CAND_COUNT,
            PEPTIDE_SIZE,
            MZ,
            MASS_THEO,
            MASS_FOUND,
            MISSED_CLEAVAGES,
            MS1COSINE_SIM,
            PPM_DIFF_MS1,
            MS1_INTENSITY,
            MONO_ISO_OFFSET,
            ISOTOPE_FOUND_CNT,
            MZ_FOUND
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
    double rescore = -1.0;

    /////
    double pVal = -1.0;
    double tTest = -1.0;
    double frameFStat = -1.0;
    double pValFrameFtest = -1.0;
    double frameError = -1.0;

    double scoreMean = -1.0;
    double scoreMedian = -1.0;
    double scoreStDev = -1.0;
    double scoreMin = -1.0;
    double scoreMax = -1.0;
    double discScoreMean = -1.0;
    double discScoreMedian = -1.0;
    double discScoreStDev = -1.0;
    double discScoreMin = -1.0;
    double discScoreMax = -1.0;
    double cosineSim = -1.0;
    double klDiv = -1.0;
    double fractionFound = -1.0;
    double mz = -1.0;
    double massTheoretical = -1.0;

    int ionsFound = -1;
    int frameCandidateCount = -1;
    int peptideSize = -1;

    double ms1CosineSim = -1.0;
    double ppmDiffMs1 = -1.0;
    double ms1Intensity = -1.0;
    double massFound = -1.0;
    double mzFound = -1.0;
    int monoOffset = -1;
    int isotopesFoundCount = -1;

    int missedCleavages= -1;

    QMap<QString, QVariant> map() override {

        using namespace PSMsReaderRowNamespace;

        return {
            {FRAME_INDEX, QVariant(frameIndex)},
            {SCAN_NUMBER, QVariant(scanNumber)},
            {CHARGE, QVariant(charge)},
            {UNIQUE_TARGET_KEY, QVariant(uniqueMsInfoScanKey)},
            {PEPTIDE_W_MODS, QVariant(peptideStringWithMods)},
            {SCORE, QVariant(score)},
            {RESCORE, QVariant(rescore)},
            {FRAME_RANK_SCORE, QVariant(frameRankScore)},
            {DISC_SCORE, QVariant(discScore)},
            {FRAME_RANK_DISC_SCORE, QVariant(frameRankDiscScore)},
            {IS_DECOY, QVariant(isDecoy)},
            {P_VAL, QVariant(pVal)},
            {T_TEST, QVariant(tTest)},
            {FRAME_F_STAT, QVariant(frameFStat)},
            {P_VAL_FRAME_F_TEST, QVariant(pValFrameFtest)},
            {FRAME_ERROR, QVariant(frameError)},
            {SCORE_MEAN , QVariant(scoreMean)},
            {SCORE_MEDIAN , QVariant(scoreMedian)},
            {SCORE_STDEV , QVariant(scoreStDev)},
            {SCORE_MAX , QVariant(scoreMin)},
            {SCORE_MAX , QVariant(scoreMax)},
            {DISC_SCORE_MEAN , QVariant(discScoreMean)},
            {DISC_SCORE_MEDIAN , QVariant(discScoreMedian)},
            {DISC_SCORE_STDEV , QVariant(discScoreStDev)},
            {DISC_SCORE_MIN , QVariant(discScoreMin)},
            {DISC_SCORE_MAX , QVariant(discScoreMax)},
            {COSINE_SIM , QVariant(cosineSim)},
            {KL_DIV , QVariant(klDiv)},
            {FRACTION_FOUND , QVariant(fractionFound)},
            {IONS_FOUND, QVariant(ionsFound)},
            {FRAME_CAND_COUNT, QVariant(frameCandidateCount)},
            {PEPTIDE_SIZE, QVariant(peptideSize)},
            {MZ, QVariant(mz)},
            {MASS_THEO, QVariant(massTheoretical)},
            {MISSED_CLEAVAGES, QVariant(missedCleavages)},
            {MS1COSINE_SIM, QVariant(ms1CosineSim)},
            {PPM_DIFF_MS1, QVariant(ppmDiffMs1)},
            {MS1_INTENSITY, QVariant(ms1Intensity)},
            {MASS_FOUND, QVariant(massFound)},
            {MONO_ISO_OFFSET, QVariant(monoOffset)},
            {ISOTOPE_FOUND_CNT, QVariant(isotopesFoundCount)},
            {MZ_FOUND, QVariant(mzFound)}

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
        rescore = dataMap.value(RESCORE).toDouble();
        frameRankScore = dataMap.value(FRAME_RANK_SCORE).toInt();
        discScore = dataMap.value(DISC_SCORE).toDouble();
        tTest = dataMap.value(T_TEST).toDouble();
        pVal = dataMap.value(P_VAL).toDouble();
        frameRankDiscScore = dataMap.value(FRAME_RANK_DISC_SCORE).toInt();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        frameFStat = dataMap.value(FRAME_F_STAT).toDouble();
        pValFrameFtest = dataMap.value(P_VAL_FRAME_F_TEST).toDouble();
        frameError = dataMap.value(FRAME_ERROR).toDouble();
        scoreMean = dataMap.value(SCORE_MEAN).toDouble();
        scoreMedian = dataMap.value(SCORE_MEDIAN).toDouble();
        scoreStDev = dataMap.value(SCORE_STDEV).toDouble();
        scoreMin = dataMap.value(SCORE_MIN).toDouble();
        scoreMax = dataMap.value(SCORE_MAX).toDouble();
        discScoreMean = dataMap.value(DISC_SCORE_MEAN).toDouble();
        discScoreMedian = dataMap.value(DISC_SCORE_MEDIAN).toDouble();
        discScoreStDev = dataMap.value(DISC_SCORE_STDEV).toDouble();
        discScoreMin = dataMap.value(DISC_SCORE_MIN).toDouble();
        discScoreMax = dataMap.value(DISC_SCORE_MAX).toDouble();
        cosineSim = dataMap.value(COSINE_SIM).toDouble();
        klDiv = dataMap.value(KL_DIV).toDouble();
        fractionFound = dataMap.value(FRACTION_FOUND).toDouble();
        ionsFound = dataMap.value(IONS_FOUND).toInt();
        frameCandidateCount = dataMap.value(FRAME_CAND_COUNT).toInt();
        peptideSize = dataMap.value(PEPTIDE_SIZE).toInt();
        missedCleavages = dataMap.value(MISSED_CLEAVAGES).toInt();
        mz = dataMap.value(MZ).toDouble();
        massTheoretical = dataMap.value(MASS_THEO).toDouble();
        ms1Intensity = dataMap.value(MS1_INTENSITY).toDouble();
        ppmDiffMs1 = dataMap.value(PPM_DIFF_MS1).toDouble();
        ms1Intensity = dataMap.value(MS1_INTENSITY).toDouble();
        massFound = dataMap.value(MASS_FOUND).toDouble();
        monoOffset = dataMap.value(MONO_ISO_OFFSET).toInt();
        isotopesFoundCount = dataMap.value(ISOTOPE_FOUND_CNT).toInt();
        mzFound = dataMap.value(MZ_FOUND).toDouble();

        ERR_RETURN
    }

};


class FILEREADERSLIB_EXPORTS PSMsReader {


public:

    PSMsReader() = default;
    ~PSMsReader() = default;

};


#endif //PYTHIADIACPP_PSMSREADER_H
