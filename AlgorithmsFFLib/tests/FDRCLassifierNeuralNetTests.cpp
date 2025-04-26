#include "FDRCLassifierNeuralNet.h"

#include "CandidateScores.h"
#include "EigenUtils.h"
#include "ParquetReader.h"

#include <QElapsedTimer>
#include <QtTest>

class FDRCLassifierNeuralNetTests : public QObject
{
    Q_OBJECT

public:
    FDRCLassifierNeuralNetTests() = default;
    ~FDRCLassifierNeuralNetTests() override = default;

private slots:

    static void playGround();


};

namespace {

    QVector<Features> getFeatures() {
        QVector<Features> featuresList = {
                CosineSimSum100,
                CosineSimSum100Top12,
                CosineSimSum100GreaterThan80,
                AllignedMaxIndexesCount,
                CosineSim100MS1,
                CosineSimSpectrumCubed,
                KlDivSpectrumCubeRoot,
                CosineSimSum45,
                CosineSimSumTop,
                CosineSimSumBottom,
                TopBottomRatio,
                TopBottomRatioNorm,
                Charge,
                ScanTimeDelta,
                ScanTimeDeltaAbs,
                ScanTimePd,
                ScanTimePdAbs,
                ScanIonCount,
                MzNorm,
                Mass,
                KlDivSpectrum,
                CosineSimSpectrum,
                CosineSim45MS1,
                CosineSim100MS1PreMono,
                CosineSim100MS1Iso1,
                CosineSim100MS1Iso2,
                PeptideLengthNorm,
                ScanTimePredicted,
                TheoFragmentCount,
                TotalIntensityLog,
                PeakShapeRatio1,
                PeakShapeRatio2,
                PeakShapeRatio3,
                ShadowsCosineSimSum,
                IRTPredicted,
                IntensityFoundMaxNorm1,
                IntensityFoundMaxNorm2,
                IntensityFoundMaxNorm3,
                IntensityFoundMaxNorm4,
                IntensityFoundMaxNorm5,
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
                AminoAcidCountK,
                AminoAcidCountL,
                AminoAcidCountM,
                AminoAcidCountN,
                AminoAcidCountP,
                AminoAcidCountQ,
                AminoAcidCountR,
                AminoAcidCountS,
                AminoAcidCountT,
                AminoAcidCountV,
                AminoAcidCountW,
                AminoAcidCountY,
                AminoAcidCountB,
                AminoAcidCountJ,
                AminoAcidCountO,
                AminoAcidCountU,
                AminoAcidCountX,
                AminoAcidCountZ,
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
                Ms1MzMeanFound100,
                Ms1MzMeanFound45,
                Ms1MzMeanFoundPreMono,
                Ms1MzMeanFoundIso1,
                Ms1MzMeanFoundIso2,
                Ms1MzMeanFound100PPM,
                Ms1MzMeanFound45PPM,
                Ms1MzMeanFoundPreMonoPPM,
                Ms1MzMeanFoundIso1PPM,
                Ms1MzMeanFoundIso2PPM,
                Ms1MzStDevFound100,
                Ms1MzStDevFound45,
                Ms1MzStDevFoundPreMono,
                Ms1MzStDevFoundIso1,
                Ms1MzStDevFoundIso2,
                Ms1IntensityFound100,
                Ms1IntensityFound45,
                Ms1IntensityFoundPreMono,
                Ms1IntensityFoundIso1,
                Ms1IntensityFoundIso2,
                CosineSimSpectrumOverTime,
                CosineSimSpectrumOverTimeCubed,
                CosineSimSpectrumStDev,
                CosineSimSum100MS1,
                // MS1Averagine,
                CosineSimSum100Window1p5X,
                CosineSimSum100Window2X,
                TotalIntensityRaw,
                TargetWindowLocation,
                DiscriminantScore,
                MatrixZeroPercentage,

                DiscScore1stRunnerUp,
                DiscScore2ndRunnerUp,
                DiscScoresCount,
                DiscScoresMean,
                DiscScoresStDev,

                // YIonSeriesMax,
                // YIonSeriesCount,
                // YIonSeriesMean,
                // YIonSeriesStd,
                // YIonSeriesTheoMax,
                // YIonSeriesTheoCount,
                // YIonSeriesTheoMean,
                // YIonSeriesTheoStd,
                // YIonSeriesMaxFoundToTheoFraction,
                // BIonSeriesMax,
                // BIonSeriesCount,
                // BIonSeriesMean,
                // BIonSeriesStd,
                // BIonSeriesTheoMax,
                // BIonSeriesTheoCount,
                // BIonSeriesTheoMean,
                // BIonSeriesTheoStd,
                // BIonSeriesMaxFoundToTheoFraction,
                // YIonSeriesCountRatio,
                // BIonSeriesCountRatio,

                CosineSimFullTheo,
                IonsFoundFractionFull,
                CosineSimFullTheoXIonsFoundFractionFull,

                Features::MzPPMMean,
                Features::FoundB,
                Features::FoundY,
                Features::FoundPercent,
                Features::MzPPMStd,
                // ImPeakLength,
                // ImAreaRatioTop2,
                // ImAreaLog10,
                // ImTheoDiff
    };
        return featuresList;
    }

}//namespace
void FDRCLassifierNeuralNetTests::playGround() {

    QSKIP("Skipping as this is for debugging");

    ERR_INIT

    const QString filePath = "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/kareNN.dat";

    QStringList proteinGroups;
    QVector<QVector<float>> xVecs;
    QVector<float> yVec;

    QFile file(filePath);
    if(file.open(QIODevice::ReadOnly)) {
        QDataStream in(&file);
        int size;
        in >> size;

        for (int i = 0; i < size; i++) {
            QVector<float> featuresArray;
            in >> featuresArray;
            featuresArray = CandidateScores::selectFeaturesArrayFeatures(
                featuresArray,
                getFeatures()
                );

            float label;
            in >> label;

            QString proteinGroup;
            in >> proteinGroup;

            xVecs.push_back(featuresArray);
            yVec.push_back(label);
            proteinGroups.push_back(proteinGroup);
        }

        file.close();
    }
    else {
        // handle error
        qDebug() << "Error reading QVector<float> from file. Reason: " << file.errorString();
    }
    qDebug() << yVec.size() << yVec.size() << proteinGroups.size() << "SIZES";

    constexpr int batchSizeMin = 500;

    const int batchSize = std::min(batchSizeMin, std::max(1, static_cast<int>(xVecs.size() / 100.0)));
    constexpr int baggingSize = 1;
    constexpr float learningRate = 0.003;
    constexpr int epochs = 3; //TODO make this settable
    FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
    e = fdrcLassifierNeuralNet.init(
            epochs,
            baggingSize,
            batchSize,
            learningRate,
            0.5,
            8
    );
    QCOMPARE(e, eNoError);

    e = fdrcLassifierNeuralNet.trainClassifier(
        xVecs,
        yVec,
        S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
        0
        );
    QCOMPARE(e, eNoError);

    QVector<float> predictions;
    e = fdrcLassifierNeuralNet.predictBaggedClassifiers(
        xVecs,
        &predictions
        );
    QCOMPARE(e, eNoError);

    struct Result {
        float y;
        float pred;
        QString proteinGroup;
    };

    QCOMPARE(yVec.size(), proteinGroups.size());
    QCOMPARE(yVec.size(), predictions.size());
    qDebug() << yVec.size() << predictions.size() << proteinGroups.size() << "SIZES";

    QVector<Result> results;
    for (int i = 0; i < yVec.size(); i++) {
        Result res;
        res.y = yVec.at(i);
        res.pred = predictions.at(i);
        res.proteinGroup = proteinGroups.at(i).isEmpty() ? "EMPY" : proteinGroups.at(i);
        results.push_back(res);
    }

    std::sort(
            results.begin(),
            results.end(),
            [](const Result &l, const Result &r){return l.pred < r.pred;}
            );

    double fdrRatio = 0.01;

    int decoyCount = 0;
    int entrapCount = 0;
    int counter = 0;
    for (int i = 0; i < results.size(); i++) {

        const Result &res = results.at(i);

        const double fdr = static_cast<double>(decoyCount) / (i + 1.0);
        if (fdr > fdrRatio) {
            counter = i;
            break;
        }

        if (res.y > 0.9) {
            decoyCount++;
            continue;
        }

        if (res.proteinGroup.contains("_ARATH") && !res.proteinGroup.contains("_HUMAN")) {
            entrapCount++;
        }

    }

    qDebug() << "fdr count" << decoyCount / static_cast<float>(counter);
    qDebug() << "entrap count" << entrapCount / static_cast<float>(counter);
    qDebug() << "PSMs found" << counter;

}

QTEST_MAIN(FDRCLassifierNeuralNetTests)

#include "FDRCLassifierNeuralNetTests.moc"
