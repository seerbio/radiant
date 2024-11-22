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

void FDRCLassifierNeuralNetTests::playGround() {

    QSKIP("Skipping as this is for debugging");

    ERR_INIT

    const QString filePath = "/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML.prqFF.candidateScores";

    // QVector<Features> featuresList = {
    //         Features::CosineSimSum100,
    //         Features::CosineSimSum100GreaterThan80,
    //         Features::AllignedMaxIndexesCount,
    //         Features::CosineSim100MS1,
    //         Features::CosineSimSpectrumCubed,
    //         Features::KlDivSpectrumCubeRoot,
    //         Features::CosineSimSum45,
    //         Features::CosineSimSumTop,
    //         Features::CosineSimSumBottom,
    //         Features::TopBottomRatio,
    //         Features::TopBottomRatioNorm,
    //         Features::Charge,
    //         Features::ScanTimeDelta,
    //         Features::ScanTimePd,
    //         Features::ScanIonCount,
    //         Features::MzNorm,
    //         Features::Mass,
    //         Features::KlDivSpectrum,
    //         Features::CosineSimSpectrum,
    //         Features::CosineSim45MS1,
    //         Features::CosineSim100MS1PreMono,
    //         Features::CosineSim100MS1Iso1,
    //         Features::CosineSim100MS1Iso2,
    //         Features::PeptideLengthNorm,
    //         Features::ScanTimePredicted,
    //         Features::TheoFragmentCount,
    //         Features::TotalIntensityLog,
    //         Features::PeakShapeRatio1,
    //         Features::PeakShapeRatio2,
    //         Features::PeakShapeRatio3,
    //         Features::ShadowsCosineSimSum,
    //         Features::IRTPredicted,
    //         Features::CosineSimToAnchor1,
    //         Features::CosineSimToAnchor2,
    //         Features::CosineSimToAnchor3,
    //         Features::CosineSimToAnchor4,
    //         Features::CosineSimToAnchor5,
    //         Features::CosineSimToAnchor6,
    //         Features::CosineSimToAnchor7,
    //         Features::CosineSimToAnchor8,
    //         Features::CosineSimToAnchor9,
    //         Features::CosineSimToAnchor10,
    //         Features::CosineSimToAnchor11,
    //         Features::CosineSimToAnchor12,
    //         Features::CosineSimShadowsToAnchor1,
    //         Features::CosineSimShadowsToAnchor2,
    //         Features::CosineSimShadowsToAnchor3,
    //         Features::CosineSimShadowsToAnchor4,
    //         Features::CosineSimShadowsToAnchor5,
    //         Features::CosineSimShadowsToAnchor6,
    //         Features::CosineSimShadowsToAnchor7,
    //         Features::CosineSimShadowsToAnchor8,
    //         Features::CosineSimShadowsToAnchor9,
    //         Features::CosineSimShadowsToAnchor10,
    //         Features::CosineSimShadowsToAnchor11,
    //         Features::CosineSimShadowsToAnchor12,
    //         Features::MzFoundMean1,
    //         Features::MzFoundMean2,
    //         Features::MzFoundMean3,
    //         Features::MzFoundMean4,
    //         Features::MzFoundMean5,
    //         Features::MzFoundMean6,
    //         Features::MzFoundMean7,
    //         Features::MzFoundMean8,
    //         Features::MzFoundMean9,
    //         Features::MzFoundMean10,
    //         Features::MzFoundMean11,
    //         Features::MzFoundMean12,
    //         Features::IntensityFoundMax1,
    //         Features::IntensityFoundMax2,
    //         Features::IntensityFoundMax3,
    //         Features::IntensityFoundMax4,
    //         Features::IntensityFoundMax5,
    //         Features::IntensityFoundMax6,
    //         Features::IntensityFoundMax7,
    //         Features::IntensityFoundMax8,
    //         Features::IntensityFoundMax9,
    //         Features::IntensityFoundMax10,
    //         Features::IntensityFoundMax11,
    //         Features::IntensityFoundMax12,
    //         Features::MzPeakLengthsNorm1,
    //         Features::MzPeakLengthsNorm2,
    //         Features::MzPeakLengthsNorm3,
    //         Features::MzPeakLengthsNorm4,
    //         Features::MzPeakLengthsNorm5,
    //         Features::MzPeakLengthsNorm6,
    //         Features::MzPeakLengthsNorm7,
    //         Features::MzPeakLengthsNorm8,
    //         Features::MzPeakLengthsNorm9,
    //         Features::MzPeakLengthsNorm10,
    //         Features::MzPeakLengthsNorm11,
    //         Features::MzPeakLengthsNorm12,
    //         Features::AminoAcidCountA,
    //         Features::AminoAcidCountC,
    //         Features::AminoAcidCountD,
    //         Features::AminoAcidCountE,
    //         Features::AminoAcidCountF,
    //         Features::AminoAcidCountG,
    //         Features::AminoAcidCountH,
    //         Features::AminoAcidCountI,
    //         Features::AminoAcidCountK,
    //         Features::AminoAcidCountL,
    //         Features::AminoAcidCountM,
    //         Features::AminoAcidCountN,
    //         Features::AminoAcidCountP,
    //         Features::AminoAcidCountQ,
    //         Features::AminoAcidCountR,
    //         Features::AminoAcidCountS,
    //         Features::AminoAcidCountT,
    //         Features::AminoAcidCountV,
    //         Features::AminoAcidCountW,
    //         Features::AminoAcidCountY,
    //         Features::AminoAcidCountB,
    //         Features::AminoAcidCountJ,
    //         Features::AminoAcidCountO,
    //         Features::AminoAcidCountU,
    //         Features::AminoAcidCountX,
    //         Features::AminoAcidCountZ,
    //         Features::MzFoundStDev1,
    //         Features::MzFoundStDev2,
    //         Features::MzFoundStDev3,
    //         Features::MzFoundStDev4,
    //         Features::MzFoundStDev5,
    //         Features::MzFoundStDev6,
    //         Features::MzFoundStDev7,
    //         Features::MzFoundStDev8,
    //         Features::MzFoundStDev9,
    //         Features::MzFoundStDev10,
    //         Features::MzFoundStDev11,
    //         Features::MzFoundStDev12,
    //         Features::AltTargetKeyIdTimeDeltaCharge1_1,
    //         Features::AltTargetKeyIdTimeDeltaCharge2_1,
    //         Features::AltTargetKeyIdTimeDeltaCharge3_1,
    //         Features::AltTargetKeyIdTimeDeltaCharge4_1,
    //         Features::Ms1MzMeanFound100,
    //         Features::Ms1MzMeanFound45,
    //         Features::Ms1MzMeanFoundPreMono,
    //         Features::Ms1MzMeanFoundIso1,
    //         Features::Ms1MzMeanFoundIso2,
    //         Features::Ms1MzMeanFound100PPM,
    //         Features::Ms1MzMeanFound45PPM,
    //         Features::Ms1MzMeanFoundPreMonoPPM,
    //         Features::Ms1MzMeanFoundIso1PPM,
    //         Features::Ms1MzMeanFoundIso2PPM,
    //         Features::Ms1MzStDevFound100,
    //         Features::Ms1MzStDevFound45,
    //         Features::Ms1MzStDevFoundPreMono,
    //         Features::Ms1MzStDevFoundIso1,
    //         Features::Ms1MzStDevFoundIso2,
    //         Features::Ms1IntensityFound100,
    //         Features::Ms1IntensityFound45,
    //         Features::Ms1IntensityFoundPreMono,
    //         Features::Ms1IntensityFoundIso1,
    //         Features::Ms1IntensityFoundIso2
    // };
    //
    // QVector<CandidateScoresReaderRow> candidateScoresReaderRows;
    // e = ParquetReader::read(
    //         filePath,
    //         &candidateScoresReaderRows
    //         );
    // QCOMPARE(e, eNoError);
    //
    // QStringList proteinGroups;
    // QVector<QVector<float>> xVecs;
    // QVector<float> yVec;
    // for (const CandidateScoresReaderRow &r : candidateScoresReaderRows) {
    //
    //     QVector<float> featuresArray = CandidateScoresReaderRow::featuresArrayFromCandidateScoresReaderRow(r);
    //     featuresArray = CandidateScores::selectFeaturesArrayFeatures(
    //             featuresArray,
    //             featuresList
    //             );
    //
    //     xVecs.push_back(featuresArray);
    //     yVec.push_back(r.isDecoy);
    //     proteinGroups.push_back(r.proteinGroup);
    // }
    //
    // Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(xVecs);
    // EigenUtils::minMaxScaleMatrix(&mat);
    //
    // const QVector<QVector<float>> xVecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);
    //
    // constexpr int baggingSize = 8;
    // const int batchSize = std::min(200, std::max(1, static_cast<int>(candidateScoresReaderRows.size() / 100.0)));
    // constexpr float learningRate = 0.003;
    // constexpr int epochs = 3; //TODO make this settable
    // FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
    // e = fdrcLassifierNeuralNet.init(
    //         epochs,
    //         baggingSize,
    //         batchSize,
    //         learningRate,
    //         8
    // );
    // QCOMPARE(e, eNoError);
    //
    // QVector<float> predictions;
    // e = fdrcLassifierNeuralNet.exec(
    //         xVecsNorm,
    //         yVec,
    //         S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
    //         0,
    //         &predictions
    // );
    // QCOMPARE(e, eNoError);
    //
    // struct Result {
    //     float y;
    //     float pred;
    //     QString proteinGroup;
    // };
    //
    // QCOMPARE(yVec.size(), proteinGroups.size());
    // QCOMPARE(yVec.size(), predictions.size());
    // qDebug() << yVec.size() << predictions.size() << proteinGroups.size() << "SIZES";
    //
    // QVector<Result> results;
    // for (int i = 0; i < yVec.size(); i++) {
    //     Result res;
    //     res.y = yVec.at(i);
    //     res.pred = predictions.at(i);
    //     res.proteinGroup = proteinGroups.at(i).isEmpty() ? "EMPY" : proteinGroups.at(i);
    //     results.push_back(res);
    // }
    //
    // std::sort(
    //         results.begin(),
    //         results.end(),
    //         [](const Result &l, const Result &r){return l.pred < r.pred;}
    //         );
    //
    // double fdrRatio = 0.01;
    //
    // int decoyCount = 0;
    // int entrapCount = 0;
    // int counter = 0;
    // for (int i = 0; i < results.size(); i++) {
    //
    //     const Result &res = results.at(i);
    //
    //     const double fdr = static_cast<double>(decoyCount) / (i + 1.0);
    //     if (fdr > fdrRatio) {
    //         counter = i;
    //         break;
    //     }
    //
    //     if (res.y > 0.9) {
    //         decoyCount++;
    //         continue;
    //     }
    //
    //     if (res.proteinGroup.contains("_ARATH") & !res.proteinGroup.contains("_HUMAN")) {
    //         entrapCount++;
    //     }
    //
    // }
    //
    // qDebug() << "fdr count" << decoyCount / static_cast<float>(counter);
    // qDebug() << "entrap count" << entrapCount / static_cast<float>(counter);
    // qDebug() << "PSMs found" << counter;

}

QTEST_MAIN(FDRCLassifierNeuralNetTests)

#include "FDRCLassifierNeuralNetTests.moc"
