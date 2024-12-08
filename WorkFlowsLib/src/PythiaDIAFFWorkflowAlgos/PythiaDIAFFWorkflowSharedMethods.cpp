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
        QVector<CandidateScores> &candidateScores,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;

    std::transform(
            candidateScores.begin(),
            candidateScores.end(),
            std::back_inserter(*candidateScoresPntrs),
            [](CandidateScores &cs){return &cs;}
            );

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
    QVector<CandidateScores> &candidateScoresVecBatch,
    const PythiaParameters &pythiaParameters,
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    QVector<CandidateScores*> *candidateScoresVecBatchPntrs,
    QMap<int, int> *fdrVsCounts,
    QVector<float> *weights
    ) {

    ERR_INIT

    const auto terminatorLogic = [&](CandidateScores *cs) {
        return cs->featuresArray[Features::CosineSimSum100] < pythiaParameters.minMs2FragCount;
    };
    const auto terminator = std::remove_if(candidateScoresVecBatchPntrs->begin(), candidateScoresVecBatchPntrs->end(), terminatorLogic);
    candidateScoresVecBatchPntrs->erase(terminator, candidateScoresVecBatchPntrs->end());

    e = buildCandidateScoresPtrs(
        candidateScoresVecBatch,
        candidateScoresVecBatchPntrs
        ); ree;

    QMap<PeptideSequenceWithModsChargeAndTargetKey , QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> peptideKeyVsTargetDecoyCandidateScoresPntrs;
    e = buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(
        *candidateScoresVecBatchPntrs,
        &peptideKeyVsTargetDecoyCandidateScoresPntrs
        ); ree;

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> featuresArrayTargetVsDecoy;
    e = DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
        peptideKeyVsTargetDecoyCandidateScoresPntrs.values().toVector(),
        useExtendedScores,
        useNeuralNetworkScores,
        &featuresArrayTargetVsDecoy
        ); ree;

    e = ErrorUtils::isEqual(
        peptideKeyVsTargetDecoyCandidateScoresPntrs.size(),
        featuresArrayTargetVsDecoy.size()
        ); ree;

    const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScorePairs
                                                        = peptideKeyVsTargetDecoyCandidateScoresPntrs.values().toVector();

    QVector<CandidateScores*> candidateScoresesPntrs;
    QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> featuresArrayTargetVsDecoyPntrs;
    e = buildProcessBatchPointers(
        candidateScorePairs,
        &featuresArrayTargetVsDecoy,
        &candidateScoresesPntrs,
        &featuresArrayTargetVsDecoyPntrs
        ); ree;

    e = DiscriminantScoretron::trainLDAClassifier(
            featuresArrayTargetVsDecoyPntrs,
            pythiaParameters.threadCount,
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

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        candidateScoresVecBatchPntrs
        ); ree;

    e = FDRCLassifierNeuralNet::outputFDRResults(
        *candidateScoresVecBatchPntrs,
        pythiaParameters.verbosity,
        fdrVsCounts
        ); ree;

    ERR_RETURN
}

namespace {

    PeptideSequenceWithModsChargeAndTargetKey buildPeptideSequenceWithModsChargeAndTargetKey(const CandidateScores* candidateScores) {
        return candidateScores->targetDecoyCandidatePair->peptideStringWithMods()
             + "|"
             + QString::number(candidateScores->targetDecoyCandidatePair->charge())
             + "|"
             + candidateScores->targetKey;
    }

}//namespace
Err PythiaDIAFFWorkflowSharedMethods::buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(
    const QVector<CandidateScores*> &candidateScores,
    QMap<PeptideSequenceWithModsChargeAndTargetKey , QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *peptideKeyVsTargetDecoyCandidateScoresPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;

    peptideKeyVsTargetDecoyCandidateScoresPntrs->clear();

    for (CandidateScores *cs: candidateScores) {

        if (cs->isDecoy) {
            (*peptideKeyVsTargetDecoyCandidateScoresPntrs)[cs->peptideSequenceWithModsChargeAndTargetKey].second = cs;
            continue;
        }
        (*peptideKeyVsTargetDecoyCandidateScoresPntrs)[cs->peptideSequenceWithModsChargeAndTargetKey].first = cs;
    }

    const bool allTargetsMatchedWithDecoy = std::all_of(
            peptideKeyVsTargetDecoyCandidateScoresPntrs->begin(),
            peptideKeyVsTargetDecoyCandidateScoresPntrs->end(),
            [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr){
                return pr.first != nullptr && pr.second != nullptr;
            }
    );
    e = ErrorUtils::isTrue(allTargetsMatchedWithDecoy); ree;

    ERR_RETURN
}

namespace {

    void filterTargetDecoyPairPointersByPrecursorMzRange(
        float precursorMzMin,
        float precursorMzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
        ) {

        const auto terminatorLogic = [precursorMzMin, precursorMzMax](const TargetDecoyCandidatePair *tdcp) {
            const float mzPrecursorTargetDecoyPair = tdcp->mz();
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

    for (const MsScanInfo &msScanInfo : msScanInfos) {

        if (msScanInfo.msLevel < 2) {
            continue;
        }

        const float precursorMzMin = msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons);
        const float precursorMzMax = msScanInfo.precursorTargetMz + (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons);

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsFiltered = targetDecoyCandidatePairs;
        filterTargetDecoyPairPointersByPrecursorMzRange(precursorMzMin, precursorMzMax, &targetDecoyCandidatePairsFiltered);

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyCandidatePairsFiltered);

        if (pythiaParameters.verbosity > 0) {
            qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyCandidatePairsFiltered.size() << "targetDecoyPairs found";
        }
    }

    ERR_RETURN
}

namespace {

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;

        // for (auto *cs : *candidateScores) {
        //     qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->targetDecoyCandidatePair->charge() << cs->isDecoy << cs->targetKey << "SDLFKJSL";
        // }

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
            row.iRTPredicted = static_cast<float>(cs->targetDecoyCandidatePair->iRt());
            row.scanTime = cs->scanTime;
            row.scanNumber = cs->scanNumber;

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
                row.mzFoundMeanVec = cs->featuresArray.mid(Features::MzFoundMean1, top6);
                row.mzFoundStDevVec = cs->featuresArray.mid(Features::MzFoundStDev1, top6);
                row.intensityFoundMaxVec = cs->featuresArray.mid(Features::IntensityFoundMax1, top6);
            }
            else {
                row.mzSearchedVec = {cs->targetDecoyCandidatePair->mz()};
                row.mzFoundMeanVec = {cs->featuresArray[Features::Ms1MzMeanFound100]};
                row.mzFoundStDevVec = {cs->featuresArray[Features::Ms1MzStDevFound100]};
                row.intensityFoundMaxVec = {cs->featuresArray[Features::Ms1IntensityFound100]};
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
