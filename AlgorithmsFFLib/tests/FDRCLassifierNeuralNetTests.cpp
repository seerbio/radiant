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

    ERR_INIT

    const QString filePath = "/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML.prqFF.candidateScores";

    QVector<CandidateScores::Features> featuresList = {
            CandidateScores::Features::CosineSimSum100,
            CandidateScores::Features::CosineSimSum100GreaterThan80,
            CandidateScores::Features::AllignedMaxIndexesCount,
            CandidateScores::Features::CosineSim100MS1,
            CandidateScores::Features::CosineSimSpectrumCubed,
            CandidateScores::Features::KlDivSpectrumCubeRoot,
            CandidateScores::Features::CosineSimSum45,
            CandidateScores::Features::CosineSimSum20,
            CandidateScores::Features::CosineSimSumTop6,
            CandidateScores::Features::CosineSimSumBottom6,
            CandidateScores::Features::TopBottomRatio,
            CandidateScores::Features::TopBottomRatioNorm,
            CandidateScores::Features::Charge,
            CandidateScores::Features::ScanTimeDelta,
            CandidateScores::Features::ScanTimeRange,
            CandidateScores::Features::ScanTimePd,
            CandidateScores::Features::ScanIonCount,
            CandidateScores::Features::MzNorm,
            CandidateScores::Features::Mass,
            CandidateScores::Features::KlDivSpectrum,
            CandidateScores::Features::CosineSimSpectrum,
            CandidateScores::Features::CosineSim45MS1,
            CandidateScores::Features::CosineSim20MS1,
            CandidateScores::Features::CosineSim100MS1PreMono,
            CandidateScores::Features::CosineSim100MS1Iso1,
            CandidateScores::Features::CosineSim100MS1Iso2,
            CandidateScores::Features::PeptideLengthNorm,
            CandidateScores::Features::ScanTimePredicted,
            CandidateScores::Features::TheoFragmentCount,
            CandidateScores::Features::TotalIntensityLog,
            CandidateScores::Features::PeakShapeRatio1,
            CandidateScores::Features::PeakShapeRatio2,
            CandidateScores::Features::PeakShapeRatio3,
            CandidateScores::Features::ShadowsCosineSimSum,
            CandidateScores::Features::IRTPredicted,
            CandidateScores::Features::CosineSimToAnchor1,
            CandidateScores::Features::CosineSimToAnchor2,
            CandidateScores::Features::CosineSimToAnchor3,
            CandidateScores::Features::CosineSimToAnchor4,
            CandidateScores::Features::CosineSimToAnchor5,
            CandidateScores::Features::CosineSimToAnchor6,
            CandidateScores::Features::CosineSimToAnchor7,
            CandidateScores::Features::CosineSimToAnchor8,
            CandidateScores::Features::CosineSimToAnchor9,
            CandidateScores::Features::CosineSimToAnchor10,
            CandidateScores::Features::CosineSimToAnchor11,
            CandidateScores::Features::CosineSimToAnchor12,
            CandidateScores::Features::CosineSimShadowsToAnchor1,
            CandidateScores::Features::CosineSimShadowsToAnchor2,
            CandidateScores::Features::CosineSimShadowsToAnchor3,
            CandidateScores::Features::CosineSimShadowsToAnchor4,
            CandidateScores::Features::CosineSimShadowsToAnchor5,
            CandidateScores::Features::CosineSimShadowsToAnchor6,
            CandidateScores::Features::CosineSimShadowsToAnchor7,
            CandidateScores::Features::CosineSimShadowsToAnchor8,
            CandidateScores::Features::CosineSimShadowsToAnchor9,
            CandidateScores::Features::CosineSimShadowsToAnchor10,
            CandidateScores::Features::CosineSimShadowsToAnchor11,
            CandidateScores::Features::CosineSimShadowsToAnchor12,
            CandidateScores::Features::ShadowsIntensityRatio1,
            CandidateScores::Features::ShadowsIntensityRatio2,
            CandidateScores::Features::ShadowsIntensityRatio3,
            CandidateScores::Features::ShadowsIntensityRatio4,
            CandidateScores::Features::ShadowsIntensityRatio5,
            CandidateScores::Features::ShadowsIntensityRatio6,
            CandidateScores::Features::ShadowsIntensityRatio7,
            CandidateScores::Features::ShadowsIntensityRatio8,
            CandidateScores::Features::ShadowsIntensityRatio9,
            CandidateScores::Features::ShadowsIntensityRatio10,
            CandidateScores::Features::ShadowsIntensityRatio11,
            CandidateScores::Features::ShadowsIntensityRatio12,
            CandidateScores::Features::MzSearched1,
            CandidateScores::Features::MzSearched2,
            CandidateScores::Features::MzSearched3,
            CandidateScores::Features::MzSearched4,
            CandidateScores::Features::MzSearched5,
            CandidateScores::Features::MzSearched6,
            CandidateScores::Features::MzSearched7,
            CandidateScores::Features::MzSearched8,
            CandidateScores::Features::MzSearched9,
            CandidateScores::Features::MzSearched10,
            CandidateScores::Features::MzSearched11,
            CandidateScores::Features::MzSearched12,
            CandidateScores::Features::TheoIntensity1,
            CandidateScores::Features::TheoIntensity2,
            CandidateScores::Features::TheoIntensity3,
            CandidateScores::Features::TheoIntensity4,
            CandidateScores::Features::TheoIntensity5,
            CandidateScores::Features::TheoIntensity6,
            CandidateScores::Features::TheoIntensity7,
            CandidateScores::Features::TheoIntensity8,
            CandidateScores::Features::TheoIntensity9,
            CandidateScores::Features::TheoIntensity10,
            CandidateScores::Features::TheoIntensity11,
            CandidateScores::Features::TheoIntensity12,
            CandidateScores::Features::MzFoundMean1,
            CandidateScores::Features::MzFoundMean2,
            CandidateScores::Features::MzFoundMean3,
            CandidateScores::Features::MzFoundMean4,
            CandidateScores::Features::MzFoundMean5,
            CandidateScores::Features::MzFoundMean6,
            CandidateScores::Features::MzFoundMean7,
            CandidateScores::Features::MzFoundMean8,
            CandidateScores::Features::MzFoundMean9,
            CandidateScores::Features::MzFoundMean10,
            CandidateScores::Features::MzFoundMean11,
            CandidateScores::Features::MzFoundMean12,
            CandidateScores::Features::IntensityFoundMax1,
            CandidateScores::Features::IntensityFoundMax2,
            CandidateScores::Features::IntensityFoundMax3,
            CandidateScores::Features::IntensityFoundMax4,
            CandidateScores::Features::IntensityFoundMax5,
            CandidateScores::Features::IntensityFoundMax6,
            CandidateScores::Features::IntensityFoundMax7,
            CandidateScores::Features::IntensityFoundMax8,
            CandidateScores::Features::IntensityFoundMax9,
            CandidateScores::Features::IntensityFoundMax10,
            CandidateScores::Features::IntensityFoundMax11,
            CandidateScores::Features::IntensityFoundMax12,
            CandidateScores::Features::MzPeakLengthsNorm1,
            CandidateScores::Features::MzPeakLengthsNorm2,
            CandidateScores::Features::MzPeakLengthsNorm3,
            CandidateScores::Features::MzPeakLengthsNorm4,
            CandidateScores::Features::MzPeakLengthsNorm5,
            CandidateScores::Features::MzPeakLengthsNorm6,
            CandidateScores::Features::MzPeakLengthsNorm7,
            CandidateScores::Features::MzPeakLengthsNorm8,
            CandidateScores::Features::MzPeakLengthsNorm9,
            CandidateScores::Features::MzPeakLengthsNorm10,
            CandidateScores::Features::MzPeakLengthsNorm11,
            CandidateScores::Features::MzPeakLengthsNorm12,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor1,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor2,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor3,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor4,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor5,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor6,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor7,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor8,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor9,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor10,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor11,
            CandidateScores::Features::ColumnApexIndexRatiosToAnchor12,
            CandidateScores::Features::AminoAcidCountA,
            CandidateScores::Features::AminoAcidCountC,
            CandidateScores::Features::AminoAcidCountD,
            CandidateScores::Features::AminoAcidCountE,
            CandidateScores::Features::AminoAcidCountF,
            CandidateScores::Features::AminoAcidCountG,
            CandidateScores::Features::AminoAcidCountH,
            CandidateScores::Features::AminoAcidCountI,
            CandidateScores::Features::AminoAcidCountK,
            CandidateScores::Features::AminoAcidCountL,
            CandidateScores::Features::AminoAcidCountM,
            CandidateScores::Features::AminoAcidCountN,
            CandidateScores::Features::AminoAcidCountP,
            CandidateScores::Features::AminoAcidCountQ,
            CandidateScores::Features::AminoAcidCountR,
            CandidateScores::Features::AminoAcidCountS,
            CandidateScores::Features::AminoAcidCountT,
            CandidateScores::Features::AminoAcidCountV,
            CandidateScores::Features::AminoAcidCountW,
            CandidateScores::Features::AminoAcidCountY,
            CandidateScores::Features::AminoAcidCountB,
            CandidateScores::Features::AminoAcidCountJ,
            CandidateScores::Features::AminoAcidCountO,
            CandidateScores::Features::AminoAcidCountU,
            CandidateScores::Features::AminoAcidCountX,
            CandidateScores::Features::AminoAcidCountZ,
            CandidateScores::Features::MzFoundStDev1,
            CandidateScores::Features::MzFoundStDev2,
            CandidateScores::Features::MzFoundStDev3,
            CandidateScores::Features::MzFoundStDev4,
            CandidateScores::Features::MzFoundStDev5,
            CandidateScores::Features::MzFoundStDev6,
            CandidateScores::Features::MzFoundStDev7,
            CandidateScores::Features::MzFoundStDev8,
            CandidateScores::Features::MzFoundStDev9,
            CandidateScores::Features::MzFoundStDev10,
            CandidateScores::Features::MzFoundStDev11,
            CandidateScores::Features::MzFoundStDev12,
            CandidateScores::Features::MzAccuracy1,
            CandidateScores::Features::MzAccuracy2,
            CandidateScores::Features::MzAccuracy3,
            CandidateScores::Features::MzAccuracy4,
            CandidateScores::Features::MzAccuracy5,
            CandidateScores::Features::MzAccuracy6,
            CandidateScores::Features::MzAccuracy7,
            CandidateScores::Features::MzAccuracy8,
            CandidateScores::Features::MzAccuracy9,
            CandidateScores::Features::MzAccuracy10,
            CandidateScores::Features::MzAccuracy11,
            CandidateScores::Features::MzAccuracy12,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_1,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2,
            CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_3,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_2,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_3,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_1,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2,
            CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3,
            CandidateScores::Features::Ms1MzMeanFound100,
            CandidateScores::Features::Ms1MzMeanFound45,
            CandidateScores::Features::Ms1MzMeanFound20,
            CandidateScores::Features::Ms1MzMeanFoundPreMono,
            CandidateScores::Features::Ms1MzMeanFoundIso1,
            CandidateScores::Features::Ms1MzMeanFoundIso2,
            CandidateScores::Features::Ms1MzMeanFound100PPM,
            CandidateScores::Features::Ms1MzMeanFound45PPM,
            CandidateScores::Features::Ms1MzMeanFound20PPM,
            CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM,
            CandidateScores::Features::Ms1MzMeanFoundIso1PPM,
            CandidateScores::Features::Ms1MzMeanFoundIso2PPM,
            CandidateScores::Features::Ms1MzStDevFound100,
            CandidateScores::Features::Ms1MzStDevFound45,
            CandidateScores::Features::Ms1MzStDevFound20,
            CandidateScores::Features::Ms1MzStDevFoundPreMono,
            CandidateScores::Features::Ms1MzStDevFoundIso1,
            CandidateScores::Features::Ms1MzStDevFoundIso2,
            CandidateScores::Features::Ms1IntensityFound100,
            CandidateScores::Features::Ms1IntensityFound45,
            CandidateScores::Features::Ms1IntensityFound20,
            CandidateScores::Features::Ms1IntensityFoundPreMono,
            CandidateScores::Features::Ms1IntensityFoundIso1,
            CandidateScores::Features::Ms1IntensityFoundIso2
    };

    QVector<CandidateScoresReaderRow> candidateScoresReaderRows;
    e = ParquetReader::read(
            filePath,
            &candidateScoresReaderRows
            );
    QCOMPARE(e, eNoError);

    QStringList proteinGroups;
    QVector<QVector<float>> xVecs;
    QVector<float> yVec;
    for (const CandidateScoresReaderRow &r : candidateScoresReaderRows) {

        QVector<float> featuresArray = CandidateScoresReaderRow::featuresArrayFromCandidateScoresReaderRow(r);
        featuresArray = CandidateScores::selectFeaturesArrayFeatures(
                featuresArray,
                featuresList
                );

        xVecs.push_back(featuresArray);
        yVec.push_back(r.isDecoy);
        proteinGroups.push_back(r.proteinGroup);
    }

    Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(xVecs);
    EigenUtils::minMaxScaleMatrix(&mat);

    const QVector<QVector<float>> xVecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

    const int baggingSize = 8;
    const int batchSize = std::min(200, std::max(1, static_cast<int>(candidateScoresReaderRows.size() / 100.0)));
    const float learningRate = 0.003;
    const int epochs = 3;
    FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
    e = fdrcLassifierNeuralNet.init(
            epochs,
            baggingSize,
            batchSize,
            learningRate
    );
    QCOMPARE(e, eNoError);

    QVector<float> predictions;
    e = fdrcLassifierNeuralNet.exec(
            xVecsNorm,
            yVec,
            S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
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

        if (res.proteinGroup.contains("_ARATH") & !res.proteinGroup.contains("_HUMAN")) {
            entrapCount++;
        }

    }

    qDebug() << "fdr count" << decoyCount / static_cast<float>(counter);
    qDebug() << "entrap count" << entrapCount / static_cast<float>(counter);
    qDebug() << "PSMs found" << counter;

}

QTEST_MAIN(FDRCLassifierNeuralNetTests)

#include "FDRCLassifierNeuralNetTests.moc"
