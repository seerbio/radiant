#include "TandemFragmentPredictotron.h"

#include <QtTest/QtTest>


class TandemFragmentPredictotronTests : public QObject
{
    Q_OBJECT

public:
    TandemFragmentPredictotronTests() = default;
    ~TandemFragmentPredictotronTests() = default;

private Q_SLOTS:
    void loadModelTest();
    void testPrediction();

};

void TandemFragmentPredictotronTests::loadModelTest()
{

    ERR_INIT

    const QString testModel
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_2.hdf5.json");

    const QString testModelFail
            = QStringLiteral("rnn_linear_charge_w_precursors_nce_22.hdf5.json");

    const int charge = 2;

    TandemFragmentPredictotron predictotron;
    e = predictotron.init(testModel, charge);
    QCOMPARE(e, eNoError);

    e = predictotron.init(testModelFail, charge);
    QCOMPARE(e, eFileError);

}

void TandemFragmentPredictotronTests::testPrediction() {

    ERR_INIT

    const int charge = 2;

    const QString testModel
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_2.hdf5.json");

    TandemFragmentPredictotron predictotron;
    e = predictotron.init(testModel, charge);
    QCOMPARE(e, eNoError);

    const QString seq = QStringLiteral("DTLMISR");
    PeptidePredictionInput ppi;
    ppi.peptideSequence = seq;
    ppi.collisionEnergy = 30;
    ppi.normalizedCollisionEnergy = 30;

    const QVector<PeptidePredictionInput> seqs = {ppi, ppi};

    QHash<QString, QVector<FragmentIon>> predictions;
    e = predictotron.batchPredictTandemSpectra(seqs, &predictions);
    QCOMPARE(e, eNoError);

    QStringList filteredResult = {
            "y1 0.28",
            "y2 0.3",
            "y3 0.59",
            "y4 1",
            "y5 0.57",
            "y4^2 0.02",
            "y5^2 0.05",
            "y1-NH3 0.04",
            "y2-NH3 0.11",
            "y3-NH3 0.04",
            "y4-NH3 0.04",
            "y3-H2O 0.02",
            "y4-H2O 0.04",
            "y5-H2O 0.02",
            "b2 0.22",
            "b3 0.03",
            "b2-H2O 0.15",
            "b3-H2O 0.22",
            "a2 0.43"
    };

    std::sort(filteredResult.begin(), filteredResult.end());

    const QString templateString = QStringLiteral("%1 %2");
    const QVector<FragmentIon> &prediction = predictions.value(predictions.keys().front());
    QCOMPARE(prediction.size(), filteredResult.size());

    int i = 0;
    for (const FragmentIon &fi : prediction) {
        qDebug() << templateString.arg(fi.ionLabel).arg(fi.intensity);
        QCOMPARE(templateString.arg(fi.ionLabel).arg(fi.intensity), filteredResult.at(i++));
    }

}


QTEST_MAIN(TandemFragmentPredictotronTests)
#include "TandemFragmentPredictotronTests.moc"
