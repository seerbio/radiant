//
// Created by Drucifer on 12/31/2021.
//

#include "Error.h"

#include "src/CommandLineParser.h"
#include "CommandLineParserUtils.h"

#include <QCoreApplication>
#include <QDirIterator>
#include <QElapsedTimer>

#include "PythiaDIAFFWorkflow.h"

using namespace Error;


namespace {

    QVector<Features> featuresCalibration() {
        QVector<Features> features = {
            CosineSimSum100,
            CosineSimSum100Top12,
            CosineSimSum100GreaterThan80,
            AllignedMaxIndexesCount,
            CosineSim100MS1, //5
            CosineSimSpectrumCubed,
            CosineSimSum45,
            CosineSimSumTop,
            CosineSimSumBottom, //10
            TopBottomRatio,
            TopBottomRatioNorm,
            ScanTimeDeltaAbs, //15
            ScanTimePdAbs,
            ScanIonCount,
            MzNorm,
            CosineSimSpectrum,
            CosineSim45MS1,//20
            CosineSim100MS1PreMono,
            CosineSim100MS1Iso1,
            CosineSim100MS1Iso2,
            PeptideLengthNorm,
            TotalIntensityLog,
            PeakShapeRatio1,
            PeakShapeRatio2,//30
            PeakShapeRatio3,
            ShadowsCosineSimSum,
            IntensityFoundMaxNorm1,
            IntensityFoundMaxNorm2,
            IntensityFoundMaxNorm3,
            IntensityFoundMaxNorm4,
            IntensityFoundMaxNorm5,//110
            IntensityFoundMaxNorm6,
            IntensityFoundMaxNorm7,
            IntensityFoundMaxNorm8,
            IntensityFoundMaxNorm9,
            IntensityFoundMaxNorm10,
            IntensityFoundMaxNorm11,
            IntensityFoundMaxNorm12,
            MzFoundStDev1,
            MzFoundStDev2,
            MzFoundStDev3,//170
            MzFoundStDev4,
            MzFoundStDev5,
            MzFoundStDev6,
            Ms1MzMeanFound100,//220
            Ms1MzMeanFound45,
            Ms1MzMeanFoundPreMono,
            Ms1MzMeanFoundIso1,
            Ms1MzMeanFoundIso2,
            Ms1MzMeanFound100PPM,
            Ms1MzMeanFound45PPM,
            Ms1MzMeanFoundPreMonoPPM,
            Ms1MzMeanFoundIso1PPM,//230
            Ms1MzMeanFoundIso2PPM,
            Ms1MzStDevFound100,
            Ms1MzStDevFound45,
            Ms1MzStDevFoundPreMono,
            Ms1MzStDevFoundIso1,
            Ms1MzStDevFoundIso2,
            CosineSimSpectrumOverTime,
            CosineSimSpectrumOverTimeCubed,
            CosineSimSpectrumStDev,
            CosineSimSum100MS1,
            MS1Averagine,
            CosineSimSum100Window1p5X,//250
            CosineSimSum100Window2X,
            TargetWindowLocationAbs,
            AlignmentIndexMean,
            AlignmentIndexStDev,
            AlignmentCombinedScore,
            MatrixZeroPercentage,
            MzPPMMeanAbs,
            MzPPMStd,
            FoundB,
            FoundY,
            FoundPercent,
            DiscScore1stRunnerUp,
            DiscScore2ndRunnerUp,
            DiscScoresCount,
            DiscScoresMean,
            DiscScoresStDev,
            YIonSeriesMax,
            YIonSeriesCount,
            YIonSeriesMean,
            YIonSeriesStd,
            YIonSeriesTheoMax,
            YIonSeriesTheoCount,
            YIonSeriesTheoMean,
            YIonSeriesTheoStd,
            YIonSeriesMaxFoundToTheoFraction,
            YIonSeriesCountRatio,
            BIonSeriesMax,
            BIonSeriesCount,
            BIonSeriesMean,
            BIonSeriesStd,
            BIonSeriesTheoMax,
            BIonSeriesTheoCount,
            BIonSeriesTheoMean,
            BIonSeriesTheoStd,
            BIonSeriesMaxFoundToTheoFraction,
            BIonSeriesCountRatio,
            CosineSimFullTheo,
            IonsFoundFractionFull,
            CosineSimFullTheoXIonsFoundFractionFull,
            CosineSimToAnchor1,
            CosineSimToAnchor2,
            CosineSimToAnchor3,
            CosineSimToAnchor4,
            CosineSimToAnchor5,
            CosineSimToAnchor6,
            CosineSimToAnchor7,
            CosineSimToAnchor8,
            CosineSimToAnchor9,
            CosineSimToAnchor10,
            CosineSimToAnchor11,
            CosineSimToAnchor12,
            MonoPreMonoRatio
        };

        return features;
    }

    QVector<Features> featuresNeuralNet() {
        QVector<Features> features = {
            CosineSimSum100,
            CosineSimSum100Top12,
            CosineSimSum100GreaterThan80,
            AllignedMaxIndexesCount,
            CosineSim100MS1, //5
            CosineSimSpectrumCubed,
            KlDivSpectrumCubeRoot,
            CosineSimSum45,
            CosineSimSumTop,
            CosineSimSumBottom, //10
            TopBottomRatio,
            TopBottomRatioNorm,
            Charge,
            ScanTimeDelta,
            ScanTimeDeltaAbs, //15
            ScanTimePd,
            ScanTimePdAbs,
            ScanIonCount,
            MzNorm,
            Mass, //20
            KlDivSpectrum,
            CosineSimSpectrum,
            CosineSim45MS1,//20
            CosineSim100MS1PreMono,
            CosineSim100MS1Iso1,
            CosineSim100MS1Iso2,
            PeptideLengthNorm,
            ScanTimePredicted,
            TheoFragmentCount,
            TotalIntensityLog,
            PeakShapeRatio1,
            PeakShapeRatio2,//30
            PeakShapeRatio3,
            ShadowsCosineSimSum,
            IRTPredicted,
            MzFoundMean1,
            MzFoundMean2,
            MzFoundMean3,
            MzFoundMean4,
            MzFoundMean5,
            MzFoundMean6,
            MzFoundMean7,//100
            MzFoundMean8,
            MzFoundMean9,
            MzFoundMean10,
            MzFoundMean11,
            MzFoundMean12,
            IntensityFoundMax1,
            IntensityFoundMax2,
            IntensityFoundMax3,
            IntensityFoundMax4,
            IntensityFoundMax5,//110
            IntensityFoundMax6,
            IntensityFoundMax7,
            IntensityFoundMax8,
            IntensityFoundMax9,
            IntensityFoundMax10,
            IntensityFoundMax11,
            IntensityFoundMax12,
            IntensityFoundMaxNorm1,
            IntensityFoundMaxNorm2,
            IntensityFoundMaxNorm3,
            IntensityFoundMaxNorm4,
            IntensityFoundMaxNorm5,//110
            IntensityFoundMaxNorm6,
            IntensityFoundMaxNorm7,
            IntensityFoundMaxNorm8,
            IntensityFoundMaxNorm9,
            IntensityFoundMaxNorm10,
            IntensityFoundMaxNorm11,
            IntensityFoundMaxNorm12,
            AminoAcidCountA,
            AminoAcidCountC,
            AminoAcidCountD,
            AminoAcidCountE,
            AminoAcidCountF,
            AminoAcidCountG,
            AminoAcidCountH,
            AminoAcidCountI,
            AminoAcidCountK,//150
            AminoAcidCountL,
            AminoAcidCountM,
            AminoAcidCountN,
            AminoAcidCountP,
            AminoAcidCountQ,
            AminoAcidCountR,
            AminoAcidCountS,
            AminoAcidCountT,
            AminoAcidCountV,
            AminoAcidCountW,//160
            AminoAcidCountY,
            AminoAcidCountB,
            AminoAcidCountJ,
            AminoAcidCountO,
            AminoAcidCountU,
            AminoAcidCountX,
            AminoAcidCountZ,
            MzFoundStDev1,
            MzFoundStDev2,
            MzFoundStDev3,//170
            MzFoundStDev4,
            MzFoundStDev5,
            MzFoundStDev6,
            AltTargetKeyIdDiscScoreChargeOG_alt,
            AltTargetKeyIdDiscScoreCharge1_1,
            AltTargetKeyIdDiscScoreCharge1_2,
            AltTargetKeyIdDiscScoreCharge2_1,
            AltTargetKeyIdDiscScoreCharge2_2,
            AltTargetKeyIdDiscScoreCharge3_1,
            AltTargetKeyIdDiscScoreCharge3_2,
            AltTargetKeyIdDiscScoreCharge4_1,
            AltTargetKeyIdDiscScoreCharge4_2,
            AltTargetKeyIdTimeDeltaCharge1_1,
            AltTargetKeyIdTimeDeltaCharge2_1,
            AltTargetKeyIdTimeDeltaCharge3_1,
            AltTargetKeyIdTimeDeltaCharge4_1,
            Ms1MzMeanFound100,//220
            Ms1MzMeanFound45,
            Ms1MzMeanFoundPreMono,
            Ms1MzMeanFoundIso1,
            Ms1MzMeanFoundIso2,
            Ms1MzMeanFound100PPM,
            Ms1MzMeanFound45PPM,
            Ms1MzMeanFoundPreMonoPPM,
            Ms1MzMeanFoundIso1PPM,//230
            Ms1MzMeanFoundIso2PPM,
            Ms1MzStDevFound100,
            Ms1MzStDevFound45,
            Ms1MzStDevFoundPreMono,
            Ms1MzStDevFoundIso1,
            Ms1MzStDevFoundIso2,
            Ms1IntensityFound100,
            Ms1IntensityFoundApex100,
            Ms1IntensityFoundApex100IM,
            Ms1IntensityFound45,
            Ms1IntensityFoundPreMono,
            Ms1IntensityFoundIso1,
            Ms1IntensityFoundIso2,
            CosineSimSpectrumOverTime,
            CosineSimSpectrumOverTimeCubed,
            CosineSimSpectrumStDev,
            CosineSimSum100MS1,
            MS1Averagine,
            CosineSimSum100Window1p5X,//250
            CosineSimSum100Window2X,
            TotalIntensityPeakHeights,
            TotalIntensityRaw,
            TargetWindowLocation,
            TargetWindowLocationAbs,
            DiscriminantScore,
            AlignmentIndexMean,
            AlignmentIndexStDev,
            AlignmentCombinedScore,
            MatrixZeroPercentage,
            MzPPMMeanAbs,
            MzPPMMean,
            MzPPMStd,
            FoundB,
            FoundY,
            FoundPercent,
            DiscScore1stRunnerUp,
            DiscScore2ndRunnerUp,
            DiscScoresCount,
            DiscScoresMean,
            DiscScoresStDev,
            YIonSeriesMax,
            YIonSeriesCount,
            YIonSeriesMean,
            YIonSeriesStd,
            YIonSeriesTheoMax,
            YIonSeriesTheoCount,
            YIonSeriesTheoMean,
            YIonSeriesTheoStd,
            YIonSeriesMaxFoundToTheoFraction,
            YIonSeriesCountRatio,
            BIonSeriesMax,
            BIonSeriesCount,
            BIonSeriesMean,
            BIonSeriesStd,
            BIonSeriesTheoMax,
            BIonSeriesTheoCount,
            BIonSeriesTheoMean,
            BIonSeriesTheoStd,
            BIonSeriesMaxFoundToTheoFraction,
            BIonSeriesCountRatio,
            CosineSimFullTheo,
            IonsFoundFractionFull,
            CosineSimFullTheoXIonsFoundFractionFull,
            CosineSimToAnchor1,
            CosineSimToAnchor2,
            CosineSimToAnchor3,
            CosineSimToAnchor4,
            CosineSimToAnchor5,
            CosineSimToAnchor6,
            CosineSimToAnchor7,//40
            CosineSimToAnchor8,
            CosineSimToAnchor9,
            CosineSimToAnchor10,
            CosineSimToAnchor11,
            CosineSimToAnchor12,
            MonoPreMonoRatio,
        };

        return features;
    }

    // Mock function to evaluate model performance
    double evaluateModel(const QSet<int>& featureSet) {
        // Dummy scoring function: Features closer to {0, 2, 4} are more important
        double score = 0;
        for (int feature : featureSet) {
            if (feature == 0 || feature == 2 || feature == 4) {
                score += 10;  // More important features contribute more
            } else {
                score += 3;   // Less important features contribute less
            }
        }
        return score;
    }

    // Sequential Backward Selection (SBS)
    QSet<int> backwardSelection(
        int num_features,
        PythiaDIAFFWorkflow *pythiaDiaFFWorkflow
        ) {

        QSet<int> selectedFeatures;
        for (int i = 0; i < num_features; ++i) {
            selectedFeatures.insert(i);
        }

        double bestScore = evaluateModel(selectedFeatures);
        QSet<int> bestSubset = selectedFeatures;

        while (selectedFeatures.size() > 1) {
            int worstFeature = -1;
            double bestNewScore = -std::numeric_limits<double>::infinity();

            for (int feature : selectedFeatures) {
                QSet<int> candidateSet = selectedFeatures;
                candidateSet.remove(feature);  // Try removing this feature

                if (const double newScore = evaluateModel(candidateSet); newScore > bestNewScore) {
                    bestNewScore = newScore;
                    worstFeature = feature;
                }
            }

            if (worstFeature == -1 || bestNewScore < bestScore) break; // Stop if no improvement
            selectedFeatures.remove(worstFeature);
            bestScore = bestNewScore;
            bestSubset = selectedFeatures;
        }

        return bestSubset;
    }

}//namespace
int main(int argc, char *argv[]) {

    ERR_INIT

    // const QString libFilePath = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib";
    // const QString parametersFilePath = "/home/andrewnichols/Desktop/Data/ConfigFiles/test_params_decoys_v2_2_1_32cpu.pythiaConfig";
    // const QVector<QString> msDataFilePaths = {
    //     "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML",
    //     // "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML"
    // };

    QElapsedTimer et;
    et.start();

    QCoreApplication app(argc, argv);
    CommandLineParser parser;

    if (!parser.validateArguments(QCoreApplication::arguments())) {
        return 1;
    }

    const CommandLineParser::CliParameters &cliParameters = parser.getCliParams();

    const QString fragLibPath = cliParameters.fragLibFilePath;
    const QString fastaFilePath = cliParameters.fastaFilePath;
    const QString pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QString msDataFile = cliParameters.msDataFile;

    int num_features = 6; // Example: 6 features

    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
            pythiaParamsFilePath,
            &pythiaParameters
            );
    if (e != eNoError) {
        qDebug() << "Error reading pythia parameters";
        return 1;
    }

    PythiaDIAFFWorkflow pythiaDiaFFWorkflow;
    e = pythiaDiaFFWorkflow.init(
            pythiaParameters,
            fragLibPath,
            fastaFilePath
    );
    if (e != eNoError) {
        qDebug() << "Error initializing Pythia Workflow Libraries";
        return 1;
    }

    QSet<int> selected = backwardSelection(
        num_features,
        &pythiaDiaFFWorkflow
        );

    qDebug() << "Selected Features: ";
    for (int feature : selected) {
        qDebug() << feature;
    }

    // e = pythiaDiaFFWorkflow.processFile(msDataFile);
    // if (e != eNoError) {
    //     qDebug() << msDataFile << "Did not run completely";
    //     return 1;
    // }
    // qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMing done in" << et.elapsed() << "mSec";

    return 0;
}
