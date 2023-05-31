//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ExtractedScansReader.h"
#include "FeatureFinderHillClusterTron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsFrameScoreVectorReader.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "TandemSpectraDeconvolvotron.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>

MsFrameScoretron::MsFrameScoretron()
: m_frameIndexMin(0)
, m_frameIndexMax(-1)
{}


namespace {

    Err setFeatureFinderParams(
            const PythiaParameters &pythiaParameters,
            FeatureFinderParameters *featureFinderParameters
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        featureFinderParameters->sigma = pythiaParameters.sigma;
        featureFinderParameters->filterLength = pythiaParameters.filterLength;
        featureFinderParameters->smoothCount = pythiaParameters.smoothCount;
        featureFinderParameters->signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        featureFinderParameters->tolerancePPM = pythiaParameters.ms2ExtractionWidthPPM;
        featureFinderParameters->skipScanCount = pythiaParameters.skipScanCount;
        featureFinderParameters->minScanCount = pythiaParameters.minScanCount;

        featureFinderParameters->printParams();

        e = ErrorUtils::isTrue(featureFinderParameters->isValid()); ree

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    params.print();

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isAboveThreshold(
            mzTargetStartStop.second,
            mzTargetStartStop.first,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;

    m_params = params;
    m_msDataFilePath = msDataFilePath;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_mzTargetStartStop = mzTargetStartStop;

    e = buildFragIonLibForTargetMz(fragLibFilePath); ree;
    e = ErrorUtils::isNotEmpty(m_fragPreds); ree;

    e = MsFrame::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop,
            &m_msFrame
    ); ree;

    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
    );ree;

    m_frameIndexMax = m_msFrame.scanCount();

    FeatureFinderParameters featureFinderParameters;
    e = setFeatureFinderParams(
            m_params,
            &featureFinderParameters
    ); ree

    e = m_featureFinderHillBuilder.init(featureFinderParameters); ree;
    e = m_featureFinderHillBuilder.buildHills(
            m_msFrame.scanNumberVsScanPoints()
            ); ree
    e = m_featureFinderHillBuilder.refineHills();

    qDebug() << "TargetKey" << uniqueMsInfoScanKey;
    qDebug() << "Scan Count" << m_msFrame.scanCount();
    qDebug() << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}

namespace {

    Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
    ){

        ERR_INIT

        const int expectedSplitSize = 2;

        const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                QString::SkipEmptyParts
        );

        e = ErrorUtils::isEqual(
                peptideSequenceChargeKeySplit.size(),
                expectedSplitSize); ree;

        *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

        e = ErrorUtils::toInt(
                peptideSequenceChargeKeySplit.back(),
                charge
                ); ree

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::buildFragIonLibForTargetMz(const QString &fragLibUri) {

    ERR_INIT

    FragLibReader fragLibReader;
    e = fragLibReader.init(fragLibUri); ree;

    const int monoOffset = 0;
    for (Charge charge = m_params.chargeStateMin; charge <= m_params.chargeStateMax; ++charge) {

        const double massStart = BiophysicalCalcs::calculateMassFromThomson(
                m_mzTargetStartStop.first,
                charge,
                monoOffset
        );

        const double massEnd = BiophysicalCalcs::calculateMassFromThomson(
                m_mzTargetStartStop.second,
                charge,
                monoOffset
        );

        QMap<PeptideSequenceChargeKey, MS2IonsSeparated> peptideSequenceChargeKeyVsMS2IonsSeparated;
        QMap<PeptideSequenceChargeKey, bool> peptideSequenceChargeKeyVsIsDecoy;
        e = fragLibReader.getMS2Ions(
                massStart,
                massEnd,
                &peptideSequenceChargeKeyVsMS2IonsSeparated,
                &peptideSequenceChargeKeyVsIsDecoy
        ); ree

        for (auto it = peptideSequenceChargeKeyVsMS2IonsSeparated.begin(); it != peptideSequenceChargeKeyVsMS2IonsSeparated.end(); it++) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            const MS2IonsSeparated &ms2IonsSeparated = it.value();

            PeptideStringWithMods peptideStringWithMods;
            Charge chargeFromPeptideSequenceChargeKey;
            e = peptideStringWithModsFromPeptideSequenceChargeKey(
                    peptideSequenceChargeKey,
                    &peptideStringWithMods,
                    &chargeFromPeptideSequenceChargeKey
            ); ree;

            if (charge != chargeFromPeptideSequenceChargeKey) {
                continue;
            }

            e = ErrorUtils::isTrue(peptideSequenceChargeKeyVsIsDecoy.contains(peptideSequenceChargeKey)); ree;

            m_fragPreds.insert(peptideStringWithMods, ms2IonsSeparated);
            m_fragPredsIsDecoy.insert(peptideStringWithMods, peptideSequenceChargeKeyVsIsDecoy.value(peptideSequenceChargeKey));
        }

    }

    ERR_RETURN
}


Err MsFrameScoretron::extractHillsForCandidtates(QString *frameHillsFilePath) {

    ERR_INIT

    const QPair<double, double> &mzTargetStartEnd = m_msFrame.precursorMzTargetStartEnd();
    qDebug() << "Processing window" << mzTargetStartEnd.first << mzTargetStartEnd.second;

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
    e = groupHillsForFrameCandidates(&pepStrWModsVsFrameIndexScoreResultOfTargets); ree

//    filterByFoundMzCount(
//            m_params.minFoundMzPeaks,
//            &m_pepStrWModsVsFrameIndexScoreResultOfTargets
//    );
//
//    *frameScoreVectorsFilePath = m_msDataFilePath + "." + m_msFrame.uniqueMsInfoScanKey() + ".frameScores";
//    e = writeFrameTargetScoreVectors(*frameScoreVectorsFilePath); ree;


    ERR_RETURN
}

namespace {




}//namespace
Err MsFrameScoretron::groupHillsForFrameCandidates(
        QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
) {

    ERR_INIT

    qDebug() << "Clustring Candidate Hills for target key:" << m_msFrame.uniqueMsInfoScanKey();
    qDebug() << "Frame size" << m_msFrame.scanCount();


    for (auto it = m_fragPreds.begin(); it != m_fragPreds.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const MS2IonsSeparated &ms2IonsTandemPred = it.value();

//        //drewholio
//        if (peptideStringWithMods != "PGQDTCQGDSGGPXTCEK") {
//            continue;
//        }
//        //drewholio

        QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> candidateFeatureFinderHills;
        e = getCandidateHills(
                ms2IonsTandemPred,
                &candidateFeatureFinderHills
                ); ree

        if (candidateFeatureFinderHills.isEmpty()){
            pepStrWModsVsFrameIndexScoreResultOfTargets->insert(peptideStringWithMods, {});
            continue;
        }

        HillsClusteringMS2 bestHillsClusteringMS2;
        e = FeatureFinderHillClusterTron::clusterHillsByFrameIndex(
                candidateFeatureFinderHills,
                m_params.mzMinDataStructure,
                m_params.mzMaxDataStructure,
                m_params.cosineSimThreshold,
                &bestHillsClusteringMS2
                ); ree;

//        qDebug() << bestHillsClusteringMS2.apexFeatureFinderHillPlus.featureFinderHill.scanNumberIndexMinMax();
//        qDebug() << bestHillsClusteringMS2.bestCosineSimScanPoints.size() << bestHillsClusteringMS2.cosineSimSum;
//        qDebug() << bestHillsClusteringMS2.bestCosineSimScanPoints;
//        break; //drewholio

    }

    ERR_RETURN
}

Err MsFrameScoretron::getCandidateHills(
        const MS2IonsSeparated &ms2IonsTandemPred,
        QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> *featureFinderHills
        ) {

    ERR_INIT


    QMap<IonIndex, QVector<FeatureFinderHill>> featureFinderHillsIonType;

    e = getHillsIonType(
            ms2IonsTandemPred.yIons,
            &featureFinderHillsIonType
            ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.bIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.y2Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y2_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.b2Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B2_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.aIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.A_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.yNH3Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_NH3_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.yH2OIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_H2O_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.bNH3Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_NH3_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.bH2OIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_H2O_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsIonType(
            ms2IonsTandemPred.precursorIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.PRECURSOR_IONS, featureFinderHillsIonType);

    ERR_RETURN
}

Err MsFrameScoretron::getHillsIonType(
        const QMap<IonIndex, MS2Ion> &ions,
        QMap<IonIndex, QVector<FeatureFinderHill>> *featureFinderHills
) {

    ERR_INIT

    for (auto it = ions.begin(); it != ions.end(); it++) {

        const IonIndex ionIndex = it.key();
        const MS2Ion &ion = it.value();

        QVector<FeatureFinderHill> ionHills;
        e = m_featureFinderHillBuilder.getHills(
                m_frameIndexMin,
                m_frameIndexMax,
                ion.mz,
                m_params.ms2ExtractionWidthPPM,
                &ionHills
        ); ree

        featureFinderHills->insert(ionIndex, ionHills);

    }

    ERR_RETURN
}
