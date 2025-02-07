//
// Created by andrewnichols on 9/28/24.
//

#include "PythiaDIAFFWorkflowSharedMethods.h"

#include "DiscriminantScoretron.h"
#include "FDRCLassifierNeuralNet.h"
#include "MsFrame.h"
#include "QValueSettertron.h"
#include "TurboXIC.h"


Err PythiaDIAFFWorkflowSharedMethods::buildUniqueMsScanInfosForProcessing(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        int numberOfUniqueScanInfosForPurpose,
        int offset,
        QVector<MsScanInfo> *uniqueMsScanInfosPurpose
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
    e = ErrorUtils::isAboveThreshold(numberOfUniqueScanInfosForPurpose, 8, ErrorUtilsParam::IncludeThreshold); ree;

    const double halfSizeScanInfos = uniqueMsScanInfos.size() / 2.0;
    const int incrementSize = std::max(1, static_cast<int>(halfSizeScanInfos / numberOfUniqueScanInfosForPurpose));

    for (int i = 0; i < uniqueMsScanInfos.size(); i += incrementSize) {

        const int distanceFromCenterRight
                = std::min(static_cast<int>(halfSizeScanInfos) + i + offset, uniqueMsScanInfos.size() - 1);

        const int distanceFromCenterLeft = std::max(static_cast<int>(halfSizeScanInfos) - i - offset, 0);

        if (i == 0) {
            uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
            continue;
        }

        uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
        uniqueMsScanInfosPurpose->push_back(uniqueMsScanInfos.at(distanceFromCenterLeft));

        if (uniqueMsScanInfosPurpose->size() >= numberOfUniqueScanInfosForPurpose) {
            break;
        }

    }

    uniqueMsScanInfosPurpose->resize(std::min(
        numberOfUniqueScanInfosForPurpose,
        uniqueMsScanInfosPurpose->size()
        ));

    ERR_RETURN
}

namespace {
    struct TurboXICLoadInput {
        QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        MzTargetKey mzTargetKey;
    };

    std::tuple<Err, MzTargetKey, TurboXIC*> buildTurboXICLogic(const TurboXICLoadInput &turboXicLoadInput) {

        ERR_INIT

        MsFrame msFrame;
        e = msFrame.init(
            turboXicLoadInput.scanNumberVsScanPoints,
            turboXicLoadInput.scanNumberVsScanTime
            ); rtee;

        auto turboXic = new TurboXIC();
        e = turboXic->init(msFrame.frameIndexVsScanPoints()); rtee;

        return {e, turboXicLoadInput.mzTargetKey, turboXic};
    }
}//namespace
Err PythiaDIAFFWorkflowSharedMethods::buildMzTargetKeyVsTurboXicPntrs(
    const QVector<MsScanInfo>& uniqueMsScanInfos,
    const QMap<ScanNumber, ScanTime>& scanNumberVsScanTime,
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
    QMap<MzTargetKey, TurboXIC*>* mzTargetKeyVsTurboXicPntr
    ) {

    ERR_INIT

            e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

    QVector<TurboXICLoadInput> turboXicLoadInputs;
    for (const MsScanInfo &msScanInfo : uniqueMsScanInfos) {

        const MzTargetKey mzTargetKey = msScanInfo.targetKey();
        e = ErrorUtils::contains(mzTargetKey, *diaTargetFrames); ree;

        TurboXICLoadInput turboXicLoadInput;
        turboXicLoadInput.scanNumberVsScanPoints = diaTargetFrames->value(mzTargetKey);
        turboXicLoadInput.scanNumberVsScanTime = scanNumberVsScanTime;
        turboXicLoadInput.mzTargetKey = mzTargetKey;

        turboXicLoadInputs.push_back(turboXicLoadInput);
    }

    QFuture<std::tuple<Err, MzTargetKey, TurboXIC*>> futures = QtConcurrent::mapped(
        turboXicLoadInputs,
        buildTurboXICLogic
        );
    futures.waitForFinished();

    for (const std::tuple<Err, MzTargetKey, TurboXIC*> &result : futures) {
        e = std::get<0>(result); ree;
        mzTargetKeyVsTurboXicPntr->insert(std::get<1>(result), std::get<2>(result));
    }

    ERR_RETURN

}

Err PythiaDIAFFWorkflowSharedMethods::buildCandidateScoresPtrs(
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScoresPairsVecBatch,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresPairsVecBatch); ree;

    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> & pr : candidateScoresPairsVecBatch) {
        candidateScoresPntrs->append({&pr.first, &pr.second});
    }

    ERR_RETURN

}

void PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(QVector<CandidateScores*> *candidateScoresPntrs) {
    std::sort(candidateScoresPntrs->rbegin(), candidateScoresPntrs->rend(), [](const CandidateScores *l, const CandidateScores *r){
        return l->discriminantScore < r->discriminantScore;
    });
}

namespace {

    Err buildProcessBatchPointers(
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScorePairs,
        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy,
        QVector<CandidateScores*> *candidateScoresesPntrs,
        QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> *featuresArrayTargetVsDecoyPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(featuresArrayTargetVsDecoy->isEmpty()); ree;

// #define CHECK_MIN_MAX_MEANS
        QVector<QVector<float>> vecs;
        for (int i = 0; i < featuresArrayTargetVsDecoy->size(); i++) {

            const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &candidatePair = candidateScorePairs.at(i);
            QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &featuresArrayPair = (*featuresArrayTargetVsDecoy)[i];
            featuresArrayTargetVsDecoyPntrs->push_back({&featuresArrayPair.first, &featuresArrayPair.second});

#ifdef CHECK_MIN_MAX_MEANS
            vecs.append({featuresArrayPair.first, featuresArrayPair.second});
#endif
            candidateScoresesPntrs->append({candidatePair.first, candidatePair.second});
        }

#ifdef CHECK_MIN_MAX_MEANS
        qDebug() << "Min max mean check";
        Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        for (int i = 0; i < mat.cols(); i++) {
            const Eigen::VectorX<float> c = mat.col(i);
            qDebug()
            << "Column" << i
            << "min" << c.minCoeff()
            << "max" << c.maxCoeff()
            << "mean" << c.mean()
            << "stdDev" << EigenUtils::calculateStDevOfVector(c)
            ;
        }
#endif

        ERR_RETURN
    }

    Err setCandidateScoresDiscriminantScore(
        const QVector<float> &discriminantScores,
        QVector<CandidateScores*> *candidateScoresesPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(discriminantScores); ree;
        e = ErrorUtils::isEqual(discriminantScores.size(), candidateScoresesPntrs->size()); ree;

        for (int i = 0; i < discriminantScores.size(); i++) {
            CandidateScores* cs = (*candidateScoresesPntrs)[i];
            cs->discriminantScore = discriminantScores.at(i);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflowSharedMethods::processBatch(
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScoresPairsVecBatch,
    const PythiaParameters &pythiaParameters,
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    QVector<CandidateScores*> *candidateScoresVecBatchPntrs,
    QMap<int, int> *fdrVsCounts,
    QVector<float> *weights
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresPairsVecBatch); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    e = buildCandidateScoresPtrs(
        candidateScoresPairsVecBatch,
        candidateScoresVecBatchPntrs
        ); ree;

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> targetDecoyCandidateScorePairsPntrs;
    e = buildTargetDecoyCandidateScorePairsPntrs(
        candidateScoresPairsVecBatch,
        &targetDecoyCandidateScorePairsPntrs
        ); ree;

    std::sort(
        targetDecoyCandidateScorePairsPntrs.begin(),
        targetDecoyCandidateScorePairsPntrs.end(),
        [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &l, const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &r) {
            return l.first->peptideSequenceWithModsChargeAndTargetKey < r.first->peptideSequenceWithModsChargeAndTargetKey;
        });

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> featuresArrayTargetVsDecoy;
    e = DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
        targetDecoyCandidateScorePairsPntrs,
        useExtendedScores,
        useNeuralNetworkScores,
        &featuresArrayTargetVsDecoy
        ); ree;

    e = ErrorUtils::isEqual(
        targetDecoyCandidateScorePairsPntrs.size(),
        featuresArrayTargetVsDecoy.size()
        ); ree;

    QVector<CandidateScores*> candidateScoresesPntrs;
    QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> featuresArrayTargetVsDecoyPntrs;
    e = buildProcessBatchPointers(
        targetDecoyCandidateScorePairsPntrs,
        &featuresArrayTargetVsDecoy,
        &candidateScoresesPntrs,
        &featuresArrayTargetVsDecoyPntrs
        ); ree;

    e = DiscriminantScoretron::trainLDAClassifier(
            featuresArrayTargetVsDecoyPntrs,
            pythiaParameters.verbosity,
            weights
            ); ree;

    QVector<FeaturesArray*> featuresArrayPntrs;
    for (const QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*> &pr : featuresArrayTargetVsDecoyPntrs) {
        featuresArrayPntrs.append({pr.first, pr.second});
    }

    QVector<float> discriminantScores;
    e = DiscriminantScoretron::applyWeights(
        *weights,
        pythiaParameters.threadCount,
        featuresArrayPntrs,
        &discriminantScores
        ); ree;

    e = setCandidateScoresDiscriminantScore(
        discriminantScores,
        &candidateScoresesPntrs
        ); ree;

    //Need to speed this up
    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &targetDecoyCandidateScorePairsPntrs
        ); ree;

    e = FDRCLassifierNeuralNet::outputFDRResults(
        *candidateScoresVecBatchPntrs,
        pythiaParameters.verbosity,
        fdrVsCounts
        ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflowSharedMethods::buildTargetDecoyCandidateScorePairsPntrs(
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &targetDecoyCandidateScorePairs,
    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidateScorePairs); ree;

    targetDecoyCandidateScorePairsPntrs->clear();

    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &pr : targetDecoyCandidateScorePairs) {
        targetDecoyCandidateScorePairsPntrs->push_back({&pr.first, &pr.second});
    }

    ERR_RETURN
}

namespace {

    void filterTargetDecoyPairPointersByPrecursorMzRange(
        float precursorMzMin,
        float precursorMzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
        ) {

        const auto terminatorLogic = [precursorMzMin, precursorMzMax](const TargetDecoyCandidatePair *tdcp) {
            const float mzPrecursorTargetDecoyPair = tdcp->mz(false);
            return !(precursorMzMin <= mzPrecursorTargetDecoyPair && mzPrecursorTargetDecoyPair <= precursorMzMax);
        };

        const auto terminator = std::remove_if(
            targetDecoyCandidatePairs->begin(),
            targetDecoyCandidatePairs->end(),
            terminatorLogic
            );

        targetDecoyCandidatePairs->erase(terminator, targetDecoyCandidatePairs->end());
    }

} //namespace
Err PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        const PythiaParameters &pythiaParameters,
        const QVector<MsScanInfo> &msScanInfos,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    mzTargetKeyVsTargetDecoyCandidatePointers->clear();

#define USE_PARALLEL_FILTERING2
#ifdef USE_PARALLEL_FILTERING2
    const auto filterBinder = std::bind(
        filterParallelLogic,
        targetDecoyCandidatePairs,
        pythiaParameters,
        std::placeholders::_1
        );

    QFuture<QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>>> futures = QtConcurrent::mapped(
        msScanInfos,
        filterBinder
        );
    futures.waitForFinished();

    for (const QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> &res: futures) {
        e = res.first; ree;
        mzTargetKeyVsTargetDecoyCandidatePointers->insert(res.second.first, res.second.second);
    }
#else
    for (const MsScanInfo &msScanInfo : msScanInfos) {
        QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> result = filterParallelLogic(
            targetDecoyCandidatePairs,
            pythiaParameters,
            msScanInfo
            );
        e = result.first; ree;
        mzTargetKeyVsTargetDecoyCandidatePointers->insert(result.second.first, result.second.second);
    }
#endif

    ERR_RETURN
}

QPair<Err, QPair<MzTargetKey ,QVector<TargetDecoyCandidatePair*>>> PythiaDIAFFWorkflowSharedMethods::filterParallelLogic(
    const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
    const PythiaParameters &pythiaParameters,
    const MsScanInfo &msScanInfo
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(msScanInfo.msLevel < 2); rree;

    const float precursorMzMin = msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons);
    const float precursorMzMax = msScanInfo.precursorTargetMz + (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons);

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsFiltered = targetDecoyCandidatePairs;
    filterTargetDecoyPairPointersByPrecursorMzRange(precursorMzMin, precursorMzMax, &targetDecoyCandidatePairsFiltered);

    if (pythiaParameters.verbosity > 0) {
        qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyCandidatePairsFiltered.size() << "targetDecoyPairs found";
    }

    return {e, {msScanInfo.targetKey(), targetDecoyCandidatePairsFiltered}};
}

namespace {

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        std::sort(
            candidateScores->rbegin(),
            candidateScores->rend(),
            [](const CandidateScores *l, const CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
            );

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;

        for (CandidateScores *cs : *candidateScores) {

            const QString key
                = cs->targetDecoyCandidatePair->peptideStringWithMods() + QString::number(cs->targetDecoyCandidatePair->charge());
            if (keyVsCandidatesFoundBest.contains(key)) {
                continue;
            }

            keyVsCandidatesFoundBest.insert(key, cs);
        }

        *candidateScores = keyVsCandidatesFoundBest.values().toVector();
    }
}
Err PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
            const MSLevelEnum &msLevel,
            const QVector<CandidateScores*> &_candidateScores,
            int verbosity,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_candidateScores); ree;

        QVector<CandidateScores*> candidateScoresFiltered = _candidateScores;
        filterDuplicateCandidateScoresByDiscriminantScore(&candidateScoresFiltered);

        if(verbosity > 0) {
            qDebug() << _candidateScores.size() << "Found for recalibartion";
            qDebug() << candidateScoresFiltered.size() << "Found for recalibartion after duplicates filtered";
        }

        constexpr int top6 = 6;
        const auto msCalibrationReaderRowsInsertLogic = [msLevel, top6](CandidateScores *cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
            row.iRTPredicted = cs->targetDecoyCandidatePair->iRt();
            row.scanTime = cs->scanTime;
            row.scanNumber = cs->scanNumber;
            row.driftTime = cs->imDriftTime;
            row.iMPredicted = cs->targetDecoyCandidatePair->iIM();

            if (msLevel == MSLevelEnum::MS2) {

                const QVector<MS2Ion> ms2Ions = cs->isDecoy
                              ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                              : cs->targetDecoyCandidatePair->ms2IonsTarget();

                QVector<float> mzSearchedVals(top6, -1.0f);
                const int maxSize = std::min(top6, ms2Ions.size());

                for (int i = 0; i < maxSize; i++) {
                    mzSearchedVals[i] = ms2Ions.at(i).mz;
                }

                row.mzSearchedVec = mzSearchedVals;
                row.mzFoundMeanVec = cs->featuresArray.mid(MzFoundMean1, top6);
                row.mzFoundStDevVec = cs->featuresArray.mid(MzFoundStDev1, top6);
                row.intensityFoundMaxVec = cs->featuresArray.mid(IntensityFoundMax1, top6);
            }
            else {
                row.mzSearchedVec = {cs->targetDecoyCandidatePair->mz(cs->isDecoy)};
                row.mzFoundMeanVec = {cs->featuresArray[Ms1MzMeanFound100]};
                row.mzFoundStDevVec = {cs->featuresArray[Ms1MzStDevFound100]};
                row.intensityFoundMaxVec = {cs->featuresArray[Ms1IntensityFound100]};
            }

            return row;
        };

        std::transform(
                candidateScoresFiltered.begin(),
                candidateScoresFiltered.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }

namespace {

    std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> buildDiaTargetFramesLoadLogic(
        const MsScanInfo &msi,
        MsReaderPointerAcc *msReaderPointerAcc,
        MsCalibratomatic *msCalibratomatic
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); rtee;

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
        e = msReaderPointerAcc->ptr->getMzTargetScanPoints(msi.targetKey(), &scanNumberVsScanPoints); rtee;

        if (msCalibratomatic->isInitCalMS2()) {
            e = msCalibratomatic->recalibrateScanPoints(
                MSLevelEnum::MS2,
                &scanNumberVsScanPoints
                ); rtee;
        }

        return {e, msi.targetKey(), scanNumberVsScanPoints};
    }

}//namespace
Err PythiaDIAFFWorkflowSharedMethods::buildDiaTargetFrames(
    const QVector<MsScanInfo> &uniqueMsScanInfos,
    MsReaderPointerAcc *msReadPointerAcc,
    MsCalibratomatic *msCalibratomatic,
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames
    ) {

    ERR_INIT

#define RUN_PARALLEL_DIATARGETFRAMES
#ifdef RUN_PARALLEL_DIATARGETFRAMES
QVector<QVector<MsScanInfo>> msScanInfosesTranced;
    e = ParallelUtils::trancheVectorForParallelization(
        uniqueMsScanInfos,
        ParallelUtils::numberOfAvailableSystemProcessors(),
        &msScanInfosesTranced
        );

    const auto trainingLogicBinder = std::bind(
        buildDiaTargetFramesLoadLogic,
        std::placeholders::_1,
        msReadPointerAcc,
        msCalibratomatic
        );

    QFuture<std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>>> futures = QtConcurrent::mapped(
        uniqueMsScanInfos,
        trainingLogicBinder
        );
    futures.waitForFinished();

    for (const std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> &res : futures) {
        e = std::get<0>(res); ree;
        (*diaTargetFrames)[std::get<1>(res)] = std::get<2>(res);
    }
#else
    for (const MsScanInfo& msi : uniqueMsScanInfos) {
        std::tuple<Err, MzTargetKey, QMap<ScanNumber, ScanPoints>> result = buildDiaTargetFramesLoadLogic(
            msi,
            msReadPointerAcc,
            msCalibratomatic
            );
        e = std::get<0>(result); ree;
        (*diaTargetFrames)[std::get<1>(result)] = std::get<2>(result);
    }
#endif

    ERR_RETURN
}
