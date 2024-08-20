#include <QtTest>
#include <QCoreApplication>

#include "ClassifierWeightsManager.h"
#include "MathUtils.h"

#include <cstdlib>
#include <boost/qvm/map_mat_vec.hpp>

class ClassifierWeightsManagerTests : public QObject {
Q_OBJECT

public:
    ClassifierWeightsManagerTests() = default;

    ~ClassifierWeightsManagerTests() override = default;

private slots:

    static void buildDataClassifier1Test();

    static void buildDataClassifier2Test();

    static void fitWeightsApplyWeightsCombinedTest();

    static void buildClassifierWeights2BigTest();


private:

    QVector<QVector<float>> m_targetsData = {
            {1.0, 2.0, 3.0, 4.0, 5.0},
            {10.0, 20.0, 30.0, 40.0, 50.0},
            {11.0, 12.0, 13.0, 14.0, 15.0}
    };

    QVector<QVector<float>> m_decoyData = {
            {10.0, 20.0, 30.0, 40.0, 50.0},
            {11.0, 12.0, 13.0, 14.0, 15.0},
            {12.0, 22.0, 32.0, 42.0, 52.0}
    };

    void buildData(
            QVector<QVector<float>*> *targets,
            QVector<QVector<float>*> *decoys
            ) {

        std::transform(
                m_targetsData.begin(),
                m_targetsData.end(),
                std::back_inserter(*targets),
                [](QVector<float> &ff){return &ff;}
                );

        std::transform(
                m_decoyData.begin(),
                m_decoyData.end(),
                std::back_inserter(*decoys),
                [](QVector<float> &ff){return &ff;}
        );
    }

};

void ClassifierWeightsManagerTests::buildDataClassifier1Test() {

    ERR_INIT

    QVector<QVector<float>*> targets;
    QVector<QVector<float>*> decoys;

    ClassifierWeightsManagerTests classifier;
    classifier.buildData(&targets, &decoys);

    QVector<QVector<float>> A;
    QVector<float> b;
    e = ClassifierWeightsManager::buildDataClassifier1(
            targets,
            decoys,
            &A,
            &b
            );
    QCOMPARE(e, eNoError);

    const QVector<QVector<float>> expectedA = {
            {21.3333, 45.3333, 69.3333, 93.3333, 117.333},
            {45.3333, 177.333, 309.333, 441.333, 573.333},
            {69.3333, 309.333, 549.333, 789.333, 1029.33},
            {93.3333, 441.333, 789.333, 1137.33, 1485.33},
            {117.333, 573.333, 1029.33, 1485.33, 1941.33}
    };
    QCOMPARE(A.size(), expectedA.size());

    const QVector<float> expectedB = {-3.66667, -6.66667, -9.66667, -12.6667, -15.6667};
    QCOMPARE(b.size(), expectedB.size());

    for (int r = 0; r < expectedA.size(); r++) {
        QCOMPARE(static_cast<int>(b.at(r)), static_cast<int>(expectedB.at(r)));
        QCOMPARE(A.at(r).size(), expectedA.at(r).size());

        for (int c = 0; c < expectedA.at(r).size(); c++) {
            QCOMPARE(static_cast<int>(A.at(r).at(c)), static_cast<int>(expectedA.at(r).at(c)));
        }
    }

}

void ClassifierWeightsManagerTests::buildDataClassifier2Test() {

    ERR_INIT

    QVector<QVector<float>*> targets;
    QVector<QVector<float>*> decoys;

    ClassifierWeightsManagerTests classifier;
    classifier.buildData(&targets, &decoys);

    QVector<QVector<float>> A;
    QVector<float> b;
    e = ClassifierWeightsManager::buildDataClassifier2(
            targets,
            decoys,
            &A,
            &b
    );
    QCOMPARE(e, eNoError);

    const QVector<QVector<float>> expectedA = {
            {15.6667, 21.6667, 27.6667, 33.6667, 39.6667},
            {21.6667, 54.6667, 87.6667, 120.667, 153.667},
            {27.6667, 87.6667, 147.667, 207.667, 267.667},
            {33.6667, 120.667, 207.667, 294.667, 381.667},
            {39.6667, 153.667, 267.667, 381.667, 495.667}
    };
    QCOMPARE(A.size(), expectedA.size());

    const QVector<float> expectedB = {-3.66667, -6.66667, -9.66667, -12.6667, -15.6667};
    QCOMPARE(b.size(), expectedB.size());

    for (int r = 0; r < expectedA.size(); r++) {
        QCOMPARE(static_cast<int>(b.at(r)), static_cast<int>(expectedB.at(r)));
        QCOMPARE(A.at(r).size(), expectedA.at(r).size());

        for (int c = 0; c < expectedA.at(r).size(); c++) {
            QCOMPARE(static_cast<int>(A.at(r).at(c)), static_cast<int>(expectedA.at(r).at(c)));
        }
    }

}

void ClassifierWeightsManagerTests::fitWeightsApplyWeightsCombinedTest() {

    ERR_INIT

    QVector<QVector<float>*> targets;
    QVector<QVector<float>*> decoys;

    ClassifierWeightsManagerTests classifier;
    classifier.buildData(&targets, &decoys);

    QVector<QVector<float>> A;
    QVector<float> b;
    e = ClassifierWeightsManager::buildDataClassifier1(
            targets,
            decoys,
            &A,
            &b
    );
    QCOMPARE(e, eNoError);

    QVector<float> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights);
    QCOMPARE(e, eNoError);

    const QVector<float> expectedWeights = {-0.190964, 0, 0, 0, 0.0034717};
    QCOMPARE(weights.size(), expectedWeights.size());

    for (int i = 0; i < weights.size(); i++) {
        QVERIFY(MathUtils::tSame(weights.at(i), expectedWeights.at(i), 3));
    }

    QVector<float> results;
    e = ClassifierWeightsManager::applyWeights(targets, weights, &results);
    QCOMPARE(e, eNoError);

    QVector<float> expectedResults = {-0.173605, -1.73605, -2.04853};
    for (int i = 0; i < results.size(); i++) {
        QVERIFY(MathUtils::tSame(results.at(i), expectedResults.at(i), 3));
    }

}

namespace {

    QVector<QVector<float>> generateRandomMatrix(
        int rows,
        int cols,
        unsigned int seed
        ) {

        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> dist(1000, 10000);


        QVector<QVector<float>> matrix(rows, QVector<float>(cols));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = static_cast<float>(dist(rng)) / 1000.0f;
            }
        }

        return matrix;
    }

}
void ClassifierWeightsManagerTests::buildClassifierWeights2BigTest() {

    ERR_INIT

    constexpr int rows = 1000000;
    constexpr int cols = 8;

    QVector<QVector<float>> matTarget  = generateRandomMatrix(rows, cols, 666);
    QVector<QVector<float>> matDecoy  = generateRandomMatrix(rows, cols, 667);

    QVector<QVector<float>*> targets;
    std::transform(
        matTarget.begin(),
        matTarget.end(),
        std::back_inserter(targets),
        [](QVector<float> &ff){return &ff;}
        );

    QVector<QVector<float>*> decoys;
    std::transform(
            matDecoy.begin(),
            matDecoy.end(),
            std::back_inserter(decoys),
            [](QVector<float> &ff){return &ff;}
    );


    QVector<QVector<float>> A;
    QVector<float> b;
    e = ClassifierWeightsManager::buildDataClassifier2(
            targets,
            decoys,
            &A,
            &b
    );
    QCOMPARE(e, eNoError);

    QVector<QVector<float>> expectedA = {
        {6.74669  , -0.00326317  , -0.00117861     , -0.01314 ,  0.000730074  ,  0.00531688  , -0.00433867  , -0.00331034},
    {-0.00326317     ,  6.75164 ,  0.000342381  ,  0.00044495  , -0.00772258  , -0.00106353  , -0.00646666  , -0.00148941},
    {-0.00117861 ,  0.000342381     ,  6.74477 , -0.000882034  , -0.00485901  , -0.00357266  ,  0.00106667  , -0.00179319},
       {-0.01314  ,  0.00044495 , -0.000882034     ,  6.75236  ,  0.00252637   ,  0.0128713  ,  0.00342081 ,  0.000237849},
    {0.000730074  , -0.00772258  , -0.00485901  ,  0.00252637     ,  6.74237  , -0.00213658  , -0.00370252  , -0.00288495},
     {0.00531688  , -0.00106353  , -0.00357266   ,  0.0128713  , -0.00213658     ,  6.75585  , -0.00204154  ,  0.00319367},
    {-0.00433867  , -0.00646666  ,  0.00106667  ,  0.00342081  , -0.00370252  , -0.00204154     ,  6.74798 ,  0.000746934},
    {-0.00331034 , -0.00148941 ,  -0.00179319 , 0.000237849 ,  -0.00288495,  0.00319367,  0.000746934,  6.74417}
    };

    for (int row = 0; row < expectedA.size(); ++row) {
        const QVector<float> &rF = A.at(row);
        const QVector<float> &rE = expectedA.at(row);

        for (int col = 0; col < rF.size(); ++col) {
            QCOMPARE(rF.at(col), rE.at(col));
        }
    }

}

QTEST_MAIN(ClassifierWeightsManagerTests)

#include "ClassifierWeightsManagerTests.moc"
