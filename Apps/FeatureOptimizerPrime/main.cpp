//
// Created by Drucifer on 12/31/2021.
//

#include "Error.h"

#include "src/CommandLineParser.h"
#include "CommandLineParserUtils.h"

#include <QCoreApplication>
#include <QDirIterator>
#include <QElapsedTimer>

#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "PythiaDIAFFWorkflow.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"

#include <random>

#include "DiscriminantScoretron.h"

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
            DiscScore1stRunnerUpDiff,
            DiscScore2ndRunnerUpDiff,
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
            DiscScore1stRunnerUpDiff,
            DiscScore2ndRunnerUpDiff,
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

    int generateRandomInt(int minInt, int maxInt) {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(minInt, maxInt - 1);

        return dist(gen);
    }

    QVector<QVector<Features>> generateRandomFeatures(
        const QVector<Features>& features,
        int featureCountMin,
        int featureCountMax,
        int numberOfFeaturesToGenerate
        ) {

        QVector<QVector<Features>> featureses;
        featureses.reserve(numberOfFeaturesToGenerate + 1);

        featureses.push_back(DiscriminantScoretron::featuresCalibration());

        const int featureCount = features.size();

        for (int i = 0; i < numberOfFeaturesToGenerate; i++) {

            const int numberOfFeatures = generateRandomInt(featureCountMin, featureCountMax);

            QVector<Features> featuresBlock(numberOfFeatures);

            QHash<int, bool> featuresHash;
            for (int i = 0; i < numberOfFeatures; i++) {
                int randomFeatureIndex = generateRandomInt(0, featureCount);
                while (featuresHash.contains(randomFeatureIndex)) {
                    randomFeatureIndex = generateRandomInt(0, featureCount);
                }

                featuresBlock[i] = features[randomFeatureIndex];
                featuresHash.insert(randomFeatureIndex, true);
            }

            featureses.push_back(featuresBlock);
        }

        return featureses;
    }

    Err optimizeCalibration(
        PythiaParameters *pythiaParameters,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
        TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
        MsReaderPointerAcc *msReaderPointerAcc,
        QVector<std::tuple<float, int, QVector<Features>>> *featuresResults
        ) {

        ERR_INIT

        constexpr int featuresCountMin = 10;
        constexpr int featuresCountMax = 30;
        constexpr int maxFeaturesToGenerate = 1e3;

        QVector<QVector<Features>> featuresCombos = generateRandomFeatures(
            featuresCalibration(),
            featuresCountMin,
            featuresCountMax,
            maxFeaturesToGenerate
            );

        int counter = 0;
        for (const QVector<Features>& features : featuresCombos) {

            QElapsedTimer et;
            et.start();

            MsCalibratomaticSettertron msCalibratomaticSettertron;
            e = msCalibratomaticSettertron.init(
                features,
                pythiaParameters,
                msReaderPointerAcc,
                targetDecoyCandidatePairManager,
                targetDecoyCandidatePairScoretron,
                false
                ); ree;

            MsCalibratomatic msCalibratomatic;
            e = msCalibratomaticSettertron.buildCalibration(&msCalibratomatic); ree;

            qDebug()
            << qPrintable(S_GLOBAL_TIMER.elapsed())
            << counter++
            << et.elapsed()
            << "mSec"
            << features
            << msCalibratomaticSettertron.batchCounter()
            << msCalibratomaticSettertron.fdrWeightedMean();

            featuresResults->push_back({
                static_cast<float>(msCalibratomaticSettertron.fdrWeightedMean()),
                msCalibratomaticSettertron.batchCounter(),
                features
            });
        }


        ERR_RETURN
    }

    Err writeResultsToFile(
        const QString &msDataFilePath,
        const QVector<std::tuple<float, int, QVector<Features>>> &featuresResults
        ) {

        ERR_INIT

        const auto longestFeatureLengthElement = std::max_element(
            featuresResults.begin(),
            featuresResults.end(),
            [](const std::tuple<float, int, QVector<Features>> &l, const std::tuple<float, int, QVector<Features>> &r) {
                // return l.second.size() < r.second.size();
                return std::get<2>(l).size() < std::get<2>(r).size();
            });

        const int longestFeatureLength = std::get<2>(*longestFeatureLengthElement).size();

        QStringList columnNames = {"Score", "Iterations"};
        for (int i = 0; i < longestFeatureLength; i++) {
            columnNames.push_back("Feature" + QString::number(i));
        }

        const QString filePathDestination = msDataFilePath + ".optResult";

        QFile file(filePathDestination);
        e = ErrorUtils::isTrue(file.open(QIODevice::WriteOnly)); ree;

        QTextStream stream(&file);

        for (const QString &s : columnNames) {
            stream << s << "\t";
        }
        stream << '\n';

        for (const std::tuple<float, int, QVector<Features>> &f : featuresResults) {
            stream << QString::number(std::get<0>(f)) << "\t";
            stream << QString::number(std::get<1>(f)) << "\t";
            for (int i : std::get<2>(f)) {
                stream << QString::number(i) << "\t";
            }
            stream << '\n';
        }

        file.close();

        ERR_RETURN
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
    const QString msDataFilePath = cliParameters.msDataFile;

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

    ///// OVERRIDES ////////////////////////////////////////
    pythiaParameters.calibrationTrainingVolume = 10000000;
    pythiaParameters.verbosity = -1;
    pythiaParameters.optimizeMode = true;
    ////////////////////////////////////////////////////////

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath);
    if (e != eNoError) {
        qDebug() << "Error reading msDataFilePath";
        return 1;
    }

    QList<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
        fragLibPath,
        pythiaParameters.useAlternativeDecoys,
        &fragLibReaderRows
        );
    if (e != eNoError) {
        qDebug() << "Error reading library";
        return 1;
    }

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            pythiaParameters,
            &fragLibReaderRows
            ); ree;
    if (e != eNoError) {
        qDebug() << "Error building target pairs";
        return 1;
    }

    TargetDecoyCandidatePairScoretron2 targetDecoyCandidatePairScoretron;
    e = targetDecoyCandidatePairScoretron.init(
        pythiaParameters,
        &msReaderPointerAcc
        ); ree;
    if (e != eNoError) {
        qDebug() << "Error building target scoretron";
        return 1;
    }

    QVector<std::tuple<float, int, QVector<Features>>> featuresResults;
    e = optimizeCalibration(
        &pythiaParameters,
        &targetDecoyCandidatePairManager,
        &targetDecoyCandidatePairScoretron,
        &msReaderPointerAcc,
        &featuresResults
        );
    if (e != eNoError) {
        qDebug() << "Error during optimization";
        return 1;
    }

    e = writeResultsToFile(msDataFilePath, featuresResults);
    if (e != eNoError) {
        qDebug() << "writing to disk failed";
        return 1;
    }



    return 0;
}
