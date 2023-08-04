//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsFrame.h"
#include "ParquetReader.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"
#include "PythiaParameterReader.h"


using namespace Error;

class TandemDeconvolverResult;
class TurboXIC;

namespace ScoredCandidateNamespace {

    const QString FREQ_PCT_SUM = QStringLiteral("frequencyPercentSum");
    const QString FREQ_PCT_SUM_POSSIBLE = QStringLiteral("frequencyPercentSumBestPossible");
    const QString KL_DIV = QStringLiteral("klDivergence");
    const QString COS_SIM_SUM = QStringLiteral("cosineSimSum");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString MASS = QStringLiteral("mass");
    const QString SCAN_NUM = QStringLiteral("scanNumber");
    const QString SCAN_TIME = QStringLiteral("scanTime");
    const QString SCAN_ION_CNT = QStringLiteral("scanIonCount");

    const QString MZ_SRCH_V = QStringLiteral("mzSearchedVec");
    const QString THEO_INTZ_V = QStringLiteral("theoIntensityVec");
    const QString MZ_FND_MEAN_V = QStringLiteral("mzFoundMeanVec");
    const QString MZ_FND_STDEV_V = QStringLiteral("mzFoundStDevVec");
    const QString INTZ_FND_MAX_V = QStringLiteral("intensityFoundMaxVec");
    const QString FRAME_IND_MAX_DIV_ANCH_V = QStringLiteral("frameIndexMaxDiffFromAnchorVec");
    const QString COS_SIM_SUM_ANCH_V = QStringLiteral("cosineSimToAnchorVec");
    const QString PK_PNT_CNT_FND_V = QStringLiteral("peakPointCountFoundVec");
    const QString FRAG_FRQ_V = QStringLiteral("fragmenFrequencyVec");
    const QString RANK_V = QStringLiteral("rankVec");

    const QStringList keysToCheck = {
            FREQ_PCT_SUM,
            FREQ_PCT_SUM_POSSIBLE,
            KL_DIV,
            COS_SIM_SUM,
            IS_DECOY,
            PEP_STR_W_MODS,
            CHARGE,
            MASS,
            SCAN_NUM,
            SCAN_TIME,
            SCAN_ION_CNT,
            MZ_SRCH_V,
            THEO_INTZ_V,
            MZ_FND_MEAN_V,
            MZ_FND_STDEV_V,
            INTZ_FND_MAX_V,
            FRAME_IND_MAX_DIV_ANCH_V,
            COS_SIM_SUM_ANCH_V,
            PK_PNT_CNT_FND_V,
            FRAG_FRQ_V,
            RANK_V
    };
}

struct FILEREADERSLIB_EXPORTS ScoredCandidate : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    double frequencyPercentSum = -1.0;
    double frequencyPercentSumBestPossible = -1.0;
    double cosineSimSum = -1.0;
    bool isDecoy = false;
    Charge charge = -1;
    double mass = -1.0;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    int scanIonCount = -1;

    QVector<int> rankVec;
    QVector<double> mzSearchedVec;
    QVector<double> theoIntensityVec;
    QVector<double> mzFoundMeanVec;
    QVector<double> mzFoundStDevVec;
    QVector<double> intensityFoundMaxVec;
    QVector<int> frameIndexMaxDiffFromAnchorVec;
    QVector<double> cosineSimToAnchorVec;
    QVector<int> peakPointCountFoundVec;
    QVector<double> fragmentFrequencyVec;

//    double discScore = -1.0;
//    double pVal = -1.0;
//    double frameError = -1.0;
//    double tTestVal = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace ScoredCandidateNamespace;

        return {
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {FREQ_PCT_SUM, QVariant(frequencyPercentSum)},
                {FREQ_PCT_SUM_POSSIBLE, QVariant(frequencyPercentSumBestPossible)},
                {COS_SIM_SUM, QVariant(cosineSimSum)},
                {IS_DECOY, QVariant(isDecoy)},
                {CHARGE, QVariant(charge)},
                {MASS, QVariant(mass)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)},
                {SCAN_ION_CNT, QVariant(scanIonCount)},
                {MZ_SRCH_V, QVariant(qVectorToQByteArray(mzSearchedVec))},
                {THEO_INTZ_V, QVariant(qVectorToQByteArray(theoIntensityVec))},
                {MZ_FND_MEAN_V, QVariant(qVectorToQByteArray(mzFoundMeanVec))},
                {MZ_FND_STDEV_V, QVariant(qVectorToQByteArray(mzFoundStDevVec))},
                {INTZ_FND_MAX_V, QVariant(qVectorToQByteArray(intensityFoundMaxVec))},
                {FRAME_IND_MAX_DIV_ANCH_V, QVariant(qVectorToQByteArray(frameIndexMaxDiffFromAnchorVec))},
                {COS_SIM_SUM_ANCH_V, QVariant(qVectorToQByteArray(cosineSimToAnchorVec))},
                {PK_PNT_CNT_FND_V, QVariant(qVectorToQByteArray(peakPointCountFoundVec))},
                {FRAG_FRQ_V, QVariant(qVectorToQByteArray(fragmentFrequencyVec))},
                {RANK_V, QVariant(qVectorToQByteArray(rankVec))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace ScoredCandidateNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frequencyPercentSum = dataMap.value(FREQ_PCT_SUM).toDouble();
        frequencyPercentSumBestPossible = dataMap.value(FREQ_PCT_SUM_POSSIBLE).toDouble();
        cosineSimSum = dataMap.value(COS_SIM_SUM).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        mass = dataMap.value(MASS).toDouble();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toDouble();
        scanIonCount = dataMap.value(SCAN_ION_CNT).toInt();

        mzSearchedVec = bytesArrayToQVector<double>(dataMap.value(MZ_SRCH_V).toByteArray());
        theoIntensityVec = bytesArrayToQVector<double>(dataMap.value(THEO_INTZ_V).toByteArray());
        mzFoundMeanVec = bytesArrayToQVector<double>(dataMap.value(MZ_FND_MEAN_V).toByteArray());
        mzFoundStDevVec = bytesArrayToQVector<double>(dataMap.value(MZ_FND_STDEV_V).toByteArray());
        intensityFoundMaxVec = bytesArrayToQVector<double>(dataMap.value(INTZ_FND_MAX_V).toByteArray());
        frameIndexMaxDiffFromAnchorVec = bytesArrayToQVector<int>(dataMap.value(FRAME_IND_MAX_DIV_ANCH_V).toByteArray());
        cosineSimToAnchorVec = bytesArrayToQVector<double>(dataMap.value(COS_SIM_SUM_ANCH_V).toByteArray());
        peakPointCountFoundVec = bytesArrayToQVector<int>(dataMap.value(PK_PNT_CNT_FND_V).toByteArray());
        fragmentFrequencyVec = bytesArrayToQVector<double>(dataMap.value(FRAG_FRQ_V).toByteArray());
        rankVec = bytesArrayToQVector<int>(dataMap.value(RANK_V).toByteArray());

        ERR_RETURN
    }
};


struct FILEREADERSLIB_EXPORTS MS2IonPeak {

    PeptideStringWithMods peptideStringWithMods;
    double mzSearched = -1.0;
    double theoIntensity = -1.0;
    int frameIndexMax = -1;
    int frameIndexStart = -1;
    int frameIndexEnd = -1;
    QVector<double> intensityVals;
    double cosineSimToAnchor = -1.0;
    int frameIndexMaxDiffFromAnchor = 0;
    int rank = -1;

    double mzFoundMean = -1.0;
    double mzFoundStDev = -1.0;
    int pointCountFound = -1;
    double fragmentFrequency = -1.0;

};


class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    friend class MissingPeptideManualTroubleshooter;
    friend class MsFrameScoretronProcessormaticTests;
    friend class MsFrameScoretronTests;

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const QString &iRTRecalibrationFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
    );

    Err scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates);


private:

    Err initMsFrame(
            const QString &msDataFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err loadFragPredsTopN(int n);

    Err buildMS2Peaks(QVector<MS2IonPeak> *ms2IonPeaks);

//    Err deconvolveScan(
//            const QVector<ScoredCandidate> &scoredCandidatesTargets,
//            const ScanPoints &scanPointsDeisotoped,
//            QVector<ScoredCandidate> *scoredCandidatesTargetsFiltered
//            );


private:

    QMap<PeptideStringWithMods, MS2IonsSeparated> m_fragPreds;
    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_fragPredsTopN;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    QMap<PeptideStringWithMods, double> m_fragPredsMass;
    QMap<PeptideStringWithMods, double> m_fragPredsIRT;
    QMap<PeptideStringWithMods, double> m_fragPredsPredictedScanTime;
    QMap<MzHashed, FrequencyPercent> m_fragmentFrequencies;

    MS2ChargeDeconvolvotron m_ms2ChargeDeconvolvotron;

    MsFrame m_msFrame;


    PythiaParameters m_params;
    QString m_msDataFilePath;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QPair<double, double> m_mzTargetStartStop;
    double m_mzStartStopMean;
};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
