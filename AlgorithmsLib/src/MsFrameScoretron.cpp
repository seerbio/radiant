//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FrameExtractsReader.h"
#include "MsFrameScoretronProcessormatic.h"
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

    QMap<PeptideStringWithMods, HillsClusteringMS2> pepStrWModsVsHillsClusteringMS2;
    e = groupHillsForFrameCandidates(&pepStrWModsVsHillsClusteringMS2); ree

    *frameHillsFilePath = m_msDataFilePath + "." + m_msFrame.uniqueMsInfoScanKey() + ".frameExtractions";
    e = writeFrameExtracts(
            pepStrWModsVsHillsClusteringMS2,
            *frameHillsFilePath
            ); ree;

    ERR_RETURN
}

Err MsFrameScoretron::groupHillsForFrameCandidates(
        QMap<PeptideStringWithMods, HillsClusteringMS2> *pepStrWModsVsHillsClusteringMS2
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

        if (bestHillsClusteringMS2.correlatedHills.size() < m_params.minFoundMzPeaks) {
            continue;
        }

        if (m_params.findIsotopologues) {
            e = findIsotopologues(&bestHillsClusteringMS2); ree
        }

        pepStrWModsVsHillsClusteringMS2->insert(peptideStringWithMods, bestHillsClusteringMS2);

    }

    ERR_RETURN
}

Err MsFrameScoretron::getCandidateHills(
        const MS2IonsSeparated &ms2IonsTandemPred,
        QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> *featureFinderHills
        ) {

    ERR_INIT

    QMap<IonIndex, QVector<FeatureFinderHill>> featureFinderHillsIonType;

    e = getHillsForIonType(
            ms2IonsTandemPred.yIons,
            &featureFinderHillsIonType
            ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.bIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.y2Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y2_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.b2Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B2_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.aIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.A_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.yNH3Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_NH3_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.yH2OIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.Y_H2O_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.bNH3Ions,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_NH3_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.bH2OIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.B_H2O_IONS, featureFinderHillsIonType);

    featureFinderHillsIonType.clear();
    e = getHillsForIonType(
            ms2IonsTandemPred.precursorIons,
            &featureFinderHillsIonType
    ); ree
    featureFinderHills->insert(S_GLOBAL_SETTINGS.PRECURSOR_IONS, featureFinderHillsIonType);

    ERR_RETURN
}

Err MsFrameScoretron::getHillsForIonType(
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

Err MsFrameScoretron::findIsotopologues(HillsClusteringMS2 *bestHillsClusteringMS2) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(bestHillsClusteringMS2->correlatedHills); ree

    for (int i = 0; i < bestHillsClusteringMS2->correlatedHills.size(); ++i) {

        FeatureFinderHillPlus &ffhPlus = bestHillsClusteringMS2->correlatedHills[i];

        double bestIsotopologueCosineSim = -1;
        Charge bestIsotopologueCharge = -1;
        double bestIsotopologueIntensity = -1.0;

        for (Charge z = 1; z <= m_params.maxIsotopologueCharge; z++) {

            const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / z;

            QVector<FeatureFinderHill> isotopologueHills;
            e = m_featureFinderHillBuilder.getHills(
                    ffhPlus.featureFinderHill.maxIntensityScanNumberIndex(),
                    ffhPlus.featureFinderHill.maxIntensityScanNumberIndex(),
                    ffhPlus.featureFinderHill.mzMean() - chargeDistance,
                    m_params.ms2ExtractionWidthPPM,
                    &isotopologueHills
            ); ree

            for (const FeatureFinderHill &ffh : isotopologueHills) {

                double cosineSimIsotopologue;
                e = FeatureFinderHillClusterTron::calculateCosineSimBetweenHills(
                        ffh,
                        ffhPlus.featureFinderHill,
                        &cosineSimIsotopologue
                        ); ree

                if (cosineSimIsotopologue > bestIsotopologueCosineSim) {
                    bestIsotopologueCosineSim = cosineSimIsotopologue;
                    bestIsotopologueCharge = z;
                    bestIsotopologueIntensity = ffh.intensityValueMax();
                }
            }
        }

        ffhPlus.bestIsotopologueCosineSim = bestIsotopologueCosineSim;
        ffhPlus.isotopologueCharge = bestIsotopologueCharge;
        ffhPlus.isotopologueIntensity = bestIsotopologueIntensity;
    }

    ERR_RETURN
}

//TOOD consider moving this monster code elsewhere, i.e., FrameExtractsReader.
namespace {

    QMap<IonType, QVector<FeatureFinderHillPlus>> separateHillsByIonType(const HillsClusteringMS2 &hillsClusteringMs2) {

        QMap<IonType, QVector<FeatureFinderHillPlus>> separatedHills;
        
        for (const FeatureFinderHillPlus &ffhp : hillsClusteringMs2.correlatedHills) {
            separatedHills[ffhp.ionType].push_back(ffhp);
        }
        
        return separatedHills;
    }

    Err unzipFeatureFinderHillPlusVector(
            const QVector<FeatureFinderHillPlus> &featureFinderHillsPlus,
            QVector<IonIndex> *ionIndexes,
            QVector<double> *mzVals,
            QVector<double> *intensities,
            QVector<double> *cosineSimToAnchors,
            QVector<int> *hillLength,
            QVector<double> *mzStd,
            QVector<double> *isotopologueCosineSim,
            QVector<int> *isotopologueCharge,
            QVector<double> *isotopologueIntensity
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featureFinderHillsPlus); ree

        for(const FeatureFinderHillPlus &ffhp : featureFinderHillsPlus) {
            ionIndexes->push_back(ffhp.ionIndex);
            mzVals->push_back(ffhp.featureFinderHill.mzMean());
            intensities->push_back(ffhp.featureFinderHill.intensityValueMax());
            cosineSimToAnchors->push_back(ffhp.cosineSimToAnchor);
            hillLength->push_back(ffhp.featureFinderHill.scanCount());
            mzStd->push_back(ffhp.featureFinderHill.mzStDev());
            isotopologueCosineSim->push_back(ffhp.bestIsotopologueCosineSim);
            isotopologueCharge->push_back(ffhp.isotopologueCharge);
            isotopologueIntensity->push_back(ffhp.isotopologueIntensity);
        }

        e = ErrorUtils::isNotEmpty(*mzVals); ree
        e = ErrorUtils::isEqual(mzVals->size(), ionIndexes->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), intensities->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), cosineSimToAnchors->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), hillLength->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), mzStd->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), isotopologueCharge->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), isotopologueCharge->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), isotopologueIntensity->size()); ree

        ERR_RETURN
    }

    Err loadFeatureFinderHillsPlusVectors(
            const HillsClusteringMS2 &hillsClusteringMs2,
            FrameExtractsReaderRow *frameExtractsReaderRow
            ) {

        ERR_INIT

        const QMap<IonType, QVector<FeatureFinderHillPlus>> separatedHills = separateHillsByIonType(hillsClusteringMs2);

        for (auto it = separatedHills.begin(); it != separatedHills.end(); it++) {

            const IonType &ionType = it.key();
            const QVector<FeatureFinderHillPlus> &featureFinderHillsPlusses = it.value();

            QVector<IonIndex> ionIndexes;
            QVector<double> mzVals;
            QVector<double> intensities;
            QVector<double> cosineSimToAnchor;
            QVector<int> hillLength;
            QVector<double> mzStd;
            QVector<double> isotopologueCosineSim;
            QVector<int> isotopologueCharge;
            QVector<double> isotopologueIntensity;
            e = unzipFeatureFinderHillPlusVector(
                    featureFinderHillsPlusses,
                    &ionIndexes,
                    &mzVals,
                    &intensities,
                    &cosineSimToAnchor,
                    &hillLength,
                    &mzStd,
                    &isotopologueCosineSim,
                    &isotopologueCharge,
                    &isotopologueIntensity
                    ); ree

            if (ionType == S_GLOBAL_SETTINGS.Y_IONS) {
                frameExtractsReaderRow->yIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->yIonMzValsActual = mzVals;
                frameExtractsReaderRow->yIonIntesitiesActual = intensities;
                frameExtractsReaderRow->yIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->yIonHillLengthActual = hillLength;
                frameExtractsReaderRow->yIonMzStdActual = mzStd;
                frameExtractsReaderRow->yIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->yIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->yIonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.Y2_IONS) {
                frameExtractsReaderRow->y2IonIndexesActual = ionIndexes;
                frameExtractsReaderRow->y2IonMzValsActual = mzVals;
                frameExtractsReaderRow->y2IonIntesitiesActual = intensities;
                frameExtractsReaderRow->y2IonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->y2IonHillLengthActual = hillLength;
                frameExtractsReaderRow->y2IonMzStdActual = mzStd;
                frameExtractsReaderRow->y2IonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->y2IonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->y2IonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.Y_NH3_IONS) {
                frameExtractsReaderRow->yNH3IonIndexesActual = ionIndexes;
                frameExtractsReaderRow->yNH3IonMzValsActual = mzVals;
                frameExtractsReaderRow->yNH3IonIntesitiesActual = intensities;
                frameExtractsReaderRow->yNH3IonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->yNH3IonHillLengthActual = hillLength;
                frameExtractsReaderRow->yNH3IonMzStdActual = mzStd;
                frameExtractsReaderRow->yNH3IonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->yNH3IonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->yNH3IonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.Y_H2O_IONS) {
                frameExtractsReaderRow->yH2OIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->yH2OIonMzValsActual = mzVals;
                frameExtractsReaderRow->yH2OIonIntesitiesActual = intensities;
                frameExtractsReaderRow->yH2OIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->yH2OIonHillLengthActual = hillLength;
                frameExtractsReaderRow->yH2OIonMzStdActual = mzStd;
                frameExtractsReaderRow->yH2OIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->yH2OIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->yH2OIonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.B_IONS) {
                frameExtractsReaderRow->bIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->bIonMzValsActual = mzVals;
                frameExtractsReaderRow->bIonIntesitiesActual = intensities;
                frameExtractsReaderRow->bIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->bIonHillLengthActual = hillLength;
                frameExtractsReaderRow->bIonMzStdActual = mzStd;
                frameExtractsReaderRow->bIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->bIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->bIonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.B2_IONS) {
                frameExtractsReaderRow->b2IonIndexesActual = ionIndexes;
                frameExtractsReaderRow->b2IonMzValsActual = mzVals;
                frameExtractsReaderRow->b2IonIntesitiesActual = intensities;
                frameExtractsReaderRow->b2IonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->b2IonHillLengthActual = hillLength;
                frameExtractsReaderRow->b2IonMzStdActual = mzStd;
                frameExtractsReaderRow->b2IonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->b2IonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->b2IonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.B_NH3_IONS) {
                frameExtractsReaderRow->bNH3IonIndexesActual = ionIndexes;
                frameExtractsReaderRow->bNH3IonMzValsActual = mzVals;
                frameExtractsReaderRow->bNH3IonIntesitiesActual = intensities;
                frameExtractsReaderRow->bNH3IonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->bNH3IonHillLengthActual = hillLength;
                frameExtractsReaderRow->bNH3IonMzStdActual = mzStd;
                frameExtractsReaderRow->bNH3IonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->bNH3IonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->bNH3IonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.B_H2O_IONS) {
                frameExtractsReaderRow->bH2OIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->bH2OIonMzValsActual = mzVals;
                frameExtractsReaderRow->bH2OIonIntesitiesActual = intensities;
                frameExtractsReaderRow->bH2OIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->bH2OIonHillLengthActual = hillLength;
                frameExtractsReaderRow->bH2OIonMzStdActual = mzStd;
                frameExtractsReaderRow->bH2OIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->bH2OIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->bH2OIonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.A_IONS) {
                frameExtractsReaderRow->aIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->aIonMzValsActual = mzVals;
                frameExtractsReaderRow->aIonIntesitiesActual = intensities;
                frameExtractsReaderRow->aIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->aIonHillLengthActual = hillLength;
                frameExtractsReaderRow->aIonMzStdActual = mzStd;
                frameExtractsReaderRow->aIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->aIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->aIonIsotopologueIntensityActual = isotopologueIntensity;
            }

            else if (ionType == S_GLOBAL_SETTINGS.PRECURSOR_IONS) {
                frameExtractsReaderRow->precursorIonIndexesActual = ionIndexes;
                frameExtractsReaderRow->precursorIonMzValsActual = mzVals;
                frameExtractsReaderRow->precursorIonIntesitiesActual = intensities;
                frameExtractsReaderRow->precursorIonHillCosineSims = cosineSimToAnchor;
                frameExtractsReaderRow->precursorIonHillLengthActual = hillLength;
                frameExtractsReaderRow->precursorIonMzStdActual = mzStd;
                frameExtractsReaderRow->precursorIonIsotopologueCosineSimActual = isotopologueCosineSim;
                frameExtractsReaderRow->precursorIonIsotopologueChargeActual = isotopologueCharge;
                frameExtractsReaderRow->precursorIonIsotopologueIntensityActual = isotopologueIntensity;
            }
        }

        ERR_RETURN
    }

    Err unzipTheoretical(
            const QMap<IonIndex, MS2Ion> &ions,
            QVector<IonIndex> *ionIndexes,
            QVector<double> *mzVals,
            QVector<double> *intensities
            ) {

        ERR_INIT

        if (ions.isEmpty()) {
            ERR_RETURN
        }

        for (auto it = ions.begin(); it != ions.end(); it++) {

            const IonIndex ionIndex = it.key();
            const MS2Ion &ms2Ion = it.value();

            ionIndexes->push_back(ionIndex);
            mzVals->push_back(ms2Ion.mz);
            intensities->push_back(ms2Ion.intensity);
        }

        e = ErrorUtils::isNotEmpty(*mzVals); ree
        e = ErrorUtils::isEqual(mzVals->size(), ionIndexes->size()); ree
        e = ErrorUtils::isEqual(mzVals->size(), intensities->size()); ree

        ERR_RETURN
    }

    Err loadTheoreticalFrags(
            const MS2IonsSeparated &ms2IonsSeparated,
            FrameExtractsReaderRow *frameExtractsReaderRow
            ){

        ERR_INIT

        e = unzipTheoretical(
                ms2IonsSeparated.yIons,
                &frameExtractsReaderRow->yIonIndexesTheo,
                &frameExtractsReaderRow->yIonMzValsTheo,
                &frameExtractsReaderRow->yIonIntesitiesTheo
                ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.y2Ions,
                &frameExtractsReaderRow->y2IonIndexesTheo,
                &frameExtractsReaderRow->y2IonMzValsTheo,
                &frameExtractsReaderRow->y2IonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.yNH3Ions,
                &frameExtractsReaderRow->yNH3IonIndexesTheo,
                &frameExtractsReaderRow->yNH3IonMzValsTheo,
                &frameExtractsReaderRow->yNH3IonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.yH2OIons,
                &frameExtractsReaderRow->yH2OIonIndexesTheo,
                &frameExtractsReaderRow->yH2OIonMzValsTheo,
                &frameExtractsReaderRow->yH2OIonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.bIons,
                &frameExtractsReaderRow->bIonIndexesTheo,
                &frameExtractsReaderRow->bIonMzValsTheo,
                &frameExtractsReaderRow->bIonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.b2Ions,
                &frameExtractsReaderRow->b2IonIndexesTheo,
                &frameExtractsReaderRow->b2IonMzValsTheo,
                &frameExtractsReaderRow->b2IonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.bNH3Ions,
                &frameExtractsReaderRow->bNH3IonIndexesTheo,
                &frameExtractsReaderRow->bNH3IonMzValsTheo,
                &frameExtractsReaderRow->bNH3IonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.bH2OIons,
                &frameExtractsReaderRow->bH2OIonIndexesTheo,
                &frameExtractsReaderRow->bH2OIonMzValsTheo,
                &frameExtractsReaderRow->bH2OIonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.aIons,
                &frameExtractsReaderRow->aIonIndexesTheo,
                &frameExtractsReaderRow->aIonMzValsTheo,
                &frameExtractsReaderRow->aIonIntesitiesTheo
        ); ree

        e = unzipTheoretical(
                ms2IonsSeparated.precursorIons,
                &frameExtractsReaderRow->precursorIonIndexesTheo,
                &frameExtractsReaderRow->precursorIonMzValsTheo,
                &frameExtractsReaderRow->precursorIonIntesitiesTheo
        ); ree

        ERR_RETURN
    }

    Err buildFrameExtractsReaderRow(
            const PeptideStringWithMods &peptideStringWithMods,
            const HillsClusteringMS2 &hillsClusteringMs2,
            const MS2IonsSeparated &ms2IonsSeparated,
            FrameExtractsReaderRow *frameExtractsReaderRow
            ) {

        ERR_INIT

        e = loadFeatureFinderHillsPlusVectors(
                hillsClusteringMs2,
                frameExtractsReaderRow
                ); ree

        e = loadTheoreticalFrags(
                ms2IonsSeparated,
                frameExtractsReaderRow
                ); ree

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::writeFrameExtracts(
        const QMap<PeptideStringWithMods, HillsClusteringMS2> &pepStrWModsVsHillsClusteringMS2,
        const QString &destinationFilePath
        ) {

    ERR_INIT

    QVector<FrameExtractsReaderRow> frameExtractsReaderRows;

    for (auto it = pepStrWModsVsHillsClusteringMS2.begin(); it != pepStrWModsVsHillsClusteringMS2.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const HillsClusteringMS2 &hillsClusteringMs2 = it.value();
        const MS2IonsSeparated &ms2IonsSeparated = m_fragPreds.value(peptideStringWithMods);

        FrameExtractsReaderRow frameExtractsReaderRow;
        frameExtractsReaderRow.isDecoy = m_fragPredsIsDecoy.value(peptideStringWithMods);
        frameExtractsReaderRow.cosineSum = hillsClusteringMs2.cosineSimSum;

        frameExtractsReaderRow.frameIndexApex
            = hillsClusteringMs2.apexFeatureFinderHillPlus.featureFinderHill.maxIntensityScanNumberIndex();

        frameExtractsReaderRow.scanNumberApex = m_msFrame.scanNumberFromFrameIndex(frameExtractsReaderRow.frameIndexApex );

        e = buildFrameExtractsReaderRow(
                peptideStringWithMods,
                hillsClusteringMs2,
                ms2IonsSeparated,
                &frameExtractsReaderRow
                ); ree

        frameExtractsReaderRows.push_back(frameExtractsReaderRow);
    }

    e = ParquetReader::write(
            frameExtractsReaderRows,
            destinationFilePath
            ); ree

    qDebug() << frameExtractsReaderRows.size() << "rows written";

    ERR_RETURN
}

