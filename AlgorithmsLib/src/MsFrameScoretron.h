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
    const QString COS_SIM = QStringLiteral("cosineSim");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString MASS = QStringLiteral("mass");
    const QString SCAN_NUM = QStringLiteral("scanNumber");
    const QString SCAN_TIME = QStringLiteral("scanTime");

    const QStringList keysToCheck = {
            FREQ_PCT_SUM,
            FREQ_PCT_SUM_POSSIBLE,
            KL_DIV,
            COS_SIM,
            IS_DECOY,
            PEP_STR_W_MODS,
            CHARGE,
            MASS,
            SCAN_NUM,
            SCAN_TIME
    };
}

struct FILEREADERSLIB_EXPORTS ScoredCandidate : public ParquetReaderInputBase {

    double frequencyPercentSum = -1.0;
    double frequencyPercentSumBestPossible = -1.0;
    double klDivergence = -1.0;
    double cosineSim = -1.0;
    bool isDecoy = false;
    PeptideStringWithMods peptideStringWithMods;
    Charge charge = -1;
    double mass = -1.0;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;

    double discScore = -1.0;
    double pVal = -1.0;
    double frameError = -1.0;
    double tTestVal = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace ScoredCandidateNamespace;

        return {
                {FREQ_PCT_SUM, QVariant(frequencyPercentSum)},
                {FREQ_PCT_SUM_POSSIBLE, QVariant(frequencyPercentSumBestPossible)},
                {KL_DIV, QVariant(klDivergence)},
                {COS_SIM, QVariant(cosineSim)},
                {IS_DECOY, QVariant(isDecoy)},
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {CHARGE, QVariant(charge)},
                {MASS, QVariant(mass)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)}
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
        klDivergence = dataMap.value(KL_DIV).toDouble();
        cosineSim = dataMap.value(COS_SIM).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        mass = dataMap.value(MASS).toDouble();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toDouble();

        ERR_RETURN
    }
};

namespace MS2IonPeakNamespace {

    const QString FRAME_IND_START = QStringLiteral("frameIndexStart");
    const QString FRAME_IND_END = QStringLiteral("frameIndexEnd");
    const QString FRAME_IND_MAX = QStringLiteral("frameIndexMax");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString MZ_SEARCHED = QStringLiteral("mzSearched");
    const QString THEO_INTZ = QStringLiteral("theoIntensity");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString COSINE_SIM_ANCH = QStringLiteral("cosineSimToAnchor");
    const QString MZ_FOUND_MEAN = QStringLiteral("mzFoundMean");
    const QString MZ_FOUND_STD = QStringLiteral("mzFoundStd");

    const QStringList keysToCheck = {
            FRAME_IND_START,
            FRAME_IND_END,
            FRAME_IND_MAX,
            PEP_STR_W_MODS,
            MZ_SEARCHED,
            THEO_INTZ,
            INTENSITY_VALS,
            COSINE_SIM_ANCH,
            MZ_FOUND_MEAN,
            MZ_FOUND_STD
    };
}

struct FILEREADERSLIB_EXPORTS MS2IonPeak : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    double mzSearched = -1.0;
    double theoIntensity = -1.0;
    int frameIndexMax = -1;
    int frameIndexStart = -1;
    int frameIndexEnd = -1;
    QVector<double> intensityVals;
    double cosineSimToAnchor = -1.0;

    double mzFoundMean = -1.0;
    double mzFoundStDev = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace MS2IonPeakNamespace;

        return {
                {FRAME_IND_START, QVariant(frameIndexStart)},
                {FRAME_IND_END, QVariant(frameIndexEnd)},
                {FRAME_IND_MAX, QVariant(frameIndexMax)},
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {MZ_SEARCHED, QVariant(mzSearched)},
                {THEO_INTZ, theoIntensity},
                {MZ_FOUND_MEAN, mzFoundMean},
                {MZ_FOUND_STD, mzFoundStDev},
                {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
                {COSINE_SIM_ANCH, QVariant(cosineSimToAnchor)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MS2IonPeakNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frameIndexStart = dataMap.value(FRAME_IND_START).toInt();
        frameIndexEnd = dataMap.value(FRAME_IND_END).toInt();
        frameIndexMax = dataMap.value(FRAME_IND_MAX).toInt();
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        mzSearched = dataMap.value(MZ_SEARCHED).toDouble();
        theoIntensity = dataMap.value(THEO_INTZ).toDouble();
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());
        cosineSimToAnchor = dataMap.value(COSINE_SIM_ANCH).toDouble();
        mzFoundMean = dataMap.value(MZ_FOUND_MEAN).toDouble();
        mzFoundStDev = dataMap.value(MZ_FOUND_STD).toDouble();

        ERR_RETURN
    }

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

    MS2ChargeDeconvolvotron m_ms2ChargeDeconvolvotron;

    MsFrame m_msFrame;


    PythiaParameters m_params;
    QString m_msDataFilePath;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QPair<double, double> m_mzTargetStartStop;
    double m_mzStartStopMean;
};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
