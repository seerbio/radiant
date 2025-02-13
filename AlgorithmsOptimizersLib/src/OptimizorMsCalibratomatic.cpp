//
// Created by andrewnichols on 2/11/25.
//

#include "OptimizorMsCalibratomatic.h"

#include "FragLibReader.h"
#include "GeneticFeatureSelection.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"


namespace {

    Err iterateMsFiles(
        const QVector<QString> &msDataFilePaths,
        const QVector<Features> &features,
        PythiaParameters *pythiaParameters,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

        ERR_INIT;

        for (const QString &msDataFilePath : msDataFilePaths) {

            MsReaderPointerAcc msReaderPointerAcc;
            e = msReaderPointerAcc.openFile(msDataFilePath); ree;

            TargetDecoyCandidatePairScoretron2 targetDecoyCandidatePairScoretron2;
            e = targetDecoyCandidatePairScoretron2.init(
                *pythiaParameters,
                &msReaderPointerAcc
                ); ree;

            MsCalibratomaticSettertron msCalibratomaticSettertron;
            e = msCalibratomaticSettertron.init(
                features,
                pythiaParameters,
                &msReaderPointerAcc,
                targetDecoyCandidatePairManager,
                &targetDecoyCandidatePairScoretron2,
                false
                ); ree;
        }

        ERR_RETURN;
    }


    QVector<Features> features() {
        const QVector<Features> baseFeatures = {
            CosineSimSum100, //0.468,
            CosineSimSum100Top12, //0.685,
            CosineSimSum100GreaterThan80, //0.846,
            AllignedMaxIndexesCount,//0.682,
            CosineSim100MS1, //0.637,
            CosineSimSpectrumCubed,//0.079,
            KlDivSpectrumCubeRoot, //0.79,
            CosineSimSum45,//0.824,
            CosineSimSumTop,//0.599,
            CosineSimSumBottom, //0.206,
            TopBottomRatio, //0.779,
            TopBottomRatioNorm,//0.438,
            ScanTimeDeltaAbs, //0.625,
            ScanTimePdAbs, //0.749,
            ScanIonCount, //0.873,
            MzNorm, // 0.382,
            KlDivSpectrum, //0.738,
            CosineSimSpectrum, //0.899,
            CosineSim45MS1, //0.843,
            CosineSim100MS1PreMono, //0.599,
            CosineSim100MS1Iso1, //0.371,
            CosineSim100MS1Iso2, //0.82,
            PeptideLengthNorm,//0.633,
            ScanTimePredicted, //0.824,
            TheoFragmentCount,//0.146,
            TotalIntensityLog, //0.404,
            PeakShapeRatio1, //0.
            PeakShapeRatio2, //0.232,
            PeakShapeRatio3, //0.165,
            ShadowsCosineSimSum, //0.049,
            IRTPredicted, //0.798,
            CosineSimSpectrumOverTime, //0.749,
            CosineSimSpectrumOverTimeCubed, //1
            CosineSimSpectrumStDev, //0.801,
            CosineSimSum100MS1,//0.206,
            MS1Averagine, //0.479,
            CosineSimSum100Window1p5X, //0.067,
            CosineSimSum100Window2X, //0.176,
            TargetWindowLocationAbs, //0.903,
            MatrixZeroPercentage,//0.659,
            MzPPMMeanAbs, //0.839,
            MzPPMStd, //0.112,
            FoundB, // 0.798,
            FoundY,//0.67,
            FoundPercent, //0.258,
            DiscScore1stRunnerUp,//0.906,
            DiscScore2ndRunnerUp, //0.745,
            DiscScoresCount,//0.206,
            DiscScoresMean,//0.285,
            DiscScoresStDev,//0.625,
            YIonSeriesMax,//0.963,
            YIonSeriesCount, //0.772,
            YIonSeriesMean, //0.401,
            YIonSeriesStd, //0.184,
            YIonSeriesTheoMax, //0.659,
            YIonSeriesTheoCount,//0.891,
            YIonSeriesTheoMean, //0.592,
            YIonSeriesTheoStd, // 0.689,
            YIonSeriesMaxFoundToTheoFraction, //0.487,
            YIonSeriesCountRatio, //0.655,
            BIonSeriesMax,//0.547,
            BIonSeriesCount, //0.052,
            BIonSeriesMean,//0.262,
            BIonSeriesStd,//0.247,
            BIonSeriesTheoMax,//0.39,
            BIonSeriesTheoCount,//0.674,
            BIonSeriesTheoMean,//0.637,
            BIonSeriesTheoStd,//0.749,
            BIonSeriesMaxFoundToTheoFraction, // 0.079,
            BIonSeriesCountRatio, //0.929,
            CosineSimFullTheo, //0.772,
            IonsFoundFractionFull, //0.682,
            CosineSimFullTheoXIonsFoundFractionFull, //0.311,
            MonoPreMonoRatio, //0.603,
        };

        return baseFeatures;
    }

}//namespace
Err OptimizorMsCalibratomatic::optimize(
    const QVector<QString> &msDataFilePaths,
    const QString &libFilePath,
    const QString &parametersFilePath,
    int maxGenerationIterations
    ) {

    ERR_INIT;

    e = ErrorUtils::isNotEmpty(msDataFilePaths); ree;
    e = ErrorUtils::isNotEmpty(libFilePath); ree;
    e = ErrorUtils::isNotEmpty(parametersFilePath); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading library";

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            libFilePath,
            &fragLibReaderRows
            ); ree;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished reading library";

    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
            parametersFilePath,
            &pythiaParameters
            );

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
        pythiaParameters,
        &fragLibReaderRows
        ); ree;

    constexpr int populationSize = 20;
    constexpr double mutationRate = 0.01;
    e = GeneticFeatureSelection::run(
        pythiaParameters,
        msDataFilePaths,
        features(),
        maxGenerationIterations,
        populationSize,
        mutationRate,
        &targetDecoyCandidatePairManager
        ); ree;

    ERR_RETURN
}
