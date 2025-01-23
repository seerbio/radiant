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
                Features::CosineSimSum100,
                Features::CosineSimSum100Top12,
                Features::CosineSimSum100GreaterThan80,
                Features::AllignedMaxIndexesCount,
                Features::CosineSim100MS1,
                Features::CosineSimSpectrumCubed,
                Features::KlDivSpectrumCubeRoot,
                Features::CosineSimSum45,
                Features::CosineSimSumTop,
                Features::CosineSimSumBottom,
                Features::TopBottomRatio,
                Features::TopBottomRatioNorm,
                Features::Charge,
                Features::ScanTimeDelta,
                Features::ScanTimeDeltaAbs,
                Features::ScanTimePd,
                Features::ScanTimePdAbs,
                Features::ScanIonCount,
                Features::MzNorm,
                Features::Mass,
                Features::KlDivSpectrum,
                Features::CosineSimSpectrum,
                Features::CosineSim45MS1,
                Features::CosineSim100MS1PreMono,
                Features::CosineSim100MS1Iso1,
                Features::CosineSim100MS1Iso2,
                Features::PeptideLengthNorm,
                Features::ScanTimePredicted,
                Features::TheoFragmentCount,
                Features::TotalIntensityLog,
                Features::PeakShapeRatio1,
                Features::PeakShapeRatio2,
                Features::PeakShapeRatio3,
                Features::ShadowsCosineSimSum,
                Features::IRTPredicted,
                Features::IntensityFoundMaxNorm1,
                Features::IntensityFoundMaxNorm2,
                Features::IntensityFoundMaxNorm3,
                Features::IntensityFoundMaxNorm4,
                Features::IntensityFoundMaxNorm5,
                Features::IntensityFoundMaxNorm6,
                Features::IntensityFoundMaxNorm7,
                Features::IntensityFoundMaxNorm8,
                Features::IntensityFoundMaxNorm9,
                Features::IntensityFoundMaxNorm10,
                Features::IntensityFoundMaxNorm11,
                Features::IntensityFoundMaxNorm12,
                Features::AminoAcidCountA,
                Features::AminoAcidCountC,
                Features::AminoAcidCountD,
                Features::AminoAcidCountE,
                Features::AminoAcidCountF,
                Features::AminoAcidCountG,
                Features::AminoAcidCountH,
                Features::AminoAcidCountI,
                Features::AminoAcidCountK,
                Features::AminoAcidCountL,
                Features::AminoAcidCountM,
                Features::AminoAcidCountN,
                Features::AminoAcidCountP,
                Features::AminoAcidCountQ,
                Features::AminoAcidCountR,
                Features::AminoAcidCountS,
                Features::AminoAcidCountT,
                Features::AminoAcidCountV,
                Features::AminoAcidCountW,
                Features::AminoAcidCountY,
                Features::AminoAcidCountB,
                Features::AminoAcidCountJ,
                Features::AminoAcidCountO,
                Features::AminoAcidCountU,
                Features::AminoAcidCountX,
                Features::AminoAcidCountZ,
                Features::AltTargetKeyIdDiscScoreChargeOG_alt,
                Features::AltTargetKeyIdDiscScoreCharge1_1,
                Features::AltTargetKeyIdDiscScoreCharge1_2,
                Features::AltTargetKeyIdDiscScoreCharge2_1,
                Features::AltTargetKeyIdDiscScoreCharge2_2,
                Features::AltTargetKeyIdDiscScoreCharge3_1,
                Features::AltTargetKeyIdDiscScoreCharge3_2,
                Features::AltTargetKeyIdDiscScoreCharge4_1,
                Features::AltTargetKeyIdDiscScoreCharge4_2,
                Features::AltTargetKeyIdTimeDeltaCharge1_1,
                Features::AltTargetKeyIdTimeDeltaCharge2_1,
                Features::AltTargetKeyIdTimeDeltaCharge3_1,
                Features::AltTargetKeyIdTimeDeltaCharge4_1,
                Features::Ms1MzMeanFound100,
                Features::Ms1MzMeanFound45,
                Features::Ms1MzMeanFoundPreMono,
                Features::Ms1MzMeanFoundIso1,
                Features::Ms1MzMeanFoundIso2,
                Features::Ms1MzMeanFound100PPM,
                Features::Ms1MzMeanFound45PPM,
                Features::Ms1MzMeanFoundPreMonoPPM,
                Features::Ms1MzMeanFoundIso1PPM,
                Features::Ms1MzMeanFoundIso2PPM,
                Features::Ms1MzStDevFound100,
                Features::Ms1MzStDevFound45,
                Features::Ms1MzStDevFoundPreMono,
                Features::Ms1MzStDevFoundIso1,
                Features::Ms1MzStDevFoundIso2,
                Features::Ms1IntensityFound100,
                Features::Ms1IntensityFound45,
                Features::Ms1IntensityFoundPreMono,
                Features::Ms1IntensityFoundIso1,
                Features::Ms1IntensityFoundIso2,
                Features::CosineSimSpectrumOverTime,
                Features::CosineSimSpectrumOverTimeCubed,
                Features::CosineSimSpectrumStDev,
                Features::CosineSimSum100MS1,
                // Features::MS1Averagine,
                Features::CosineSimSum100Window1p5X,
                Features::CosineSimSum100Window2X,
                Features::TotalIntensityRaw,
                Features::TargetWindowLocation,
                Features::DiscriminantScore,
                Features::MatrixZeroPercentage,

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

    // QSKIP("Skipping as this is for debugging");

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
