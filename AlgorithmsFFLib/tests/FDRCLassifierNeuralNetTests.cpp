#include "FDRCLassifierNeuralNet.h"

#include "CandidateScores.h"
#include "EigenUtils.h"
#include "ParallelUtils.h"
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
        const QVector<float> featuresArray = CandidateScoresReaderRow::featuresArrayFromCandidateScoresReaderRow(r);
        xVecs.push_back(featuresArray);
        yVec.push_back(r.isDecoy);
        proteinGroups.push_back(r.proteinGroup);
        qDebug() << r.proteinGroup;
    }

    Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(xVecs);
    EigenUtils::minMaxScaleMatrix(&mat);

    const QVector<QVector<float>> xVecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

    const int baggingSize = 12;
    const int batchSize = std::min(200, std::max(1, static_cast<int>(candidateScoresReaderRows.size() / 100.0)));
    const float learningRate = 0.003;
    const int epochs = 1;
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

    QVector<QPair<float, float>> results;
    e = ParallelUtils::zip(
            yVec,
            predictions,
            &results
            );
    QCOMPARE(e, eNoError);

}

QTEST_MAIN(FDRCLassifierNeuralNetTests)

#include "FDRCLassifierNeuralNetTests.moc"
