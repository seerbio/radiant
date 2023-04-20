#include "NearestNeighborsSearch.h"

#include <QElapsedTimer>
#include <QtTest/QtTest>


class NearestNeighborsSearchTests : public QObject
{
    Q_OBJECT

public:
    NearestNeighborsSearchTests() = default;
    ~NearestNeighborsSearchTests() override = default;

private Q_SLOTS:
    void initTest();
    void kNearestNeighborsSearchTest();
    void radiusSearchTest();
    void readWriteNearestNeighborsTest();


private:

    QVector<QPair<double, QVector<double>>> buildTestPoints() {

        QVector<QPair<double, QVector<double>>> tp;

        for (int i = 0; i < 10; i++) {
            double p = 1.0 + i;
            tp.push_back({p, {p, p}});
        }

        return tp;
    }

};


void NearestNeighborsSearchTests::initTest()
{

    ERR_INIT

    const QVector<QPair<double, QVector<double>>> testPoints = buildTestPoints();

    NearestNeighborsSearch seeker;
    e = seeker.init({});
    QCOMPARE(e, eError);

    e = seeker.init(testPoints);
    QCOMPARE(e, eNoError);

}

void NearestNeighborsSearchTests::kNearestNeighborsSearchTest() {

    ERR_INIT

    const QVector<QPair<double, QVector<double>>> testPoints = buildTestPoints();

    NearestNeighborsSearch seeker;
    e = seeker.init(testPoints);
    QCOMPARE(e, eNoError);

    const QVector<QVector<double>> searchPoint = {{5.0, 5.0}, {5.0, 5.0}};

    QVector<NNSearchResult> result;
    e = seeker.kNearestNeighborsSearch(
            searchPoint,
            3,
            &result
            );
    QCOMPARE(e, eNoError);

    const int expectedResultSize = 2;
    QCOMPARE(result.size(), expectedResultSize);

    std::vector<long> expectedIndexes = {4, 3, 5};
    std::vector<double> expectedValues = {5, 4, 6};
    std::vector<double> expectedDistancesSquared = {0, 2, 2};

    const NNSearchResult &res = result.front();
    QCOMPARE(res.indexes, expectedIndexes);
    QCOMPARE(res.values, expectedValues);
    QCOMPARE(res.distancesSquared, expectedDistancesSquared);
}

void NearestNeighborsSearchTests::radiusSearchTest() {

    ERR_INIT

    const QVector<QPair<double, QVector<double>>> testPoints = buildTestPoints();

    NearestNeighborsSearch seeker;
    e = seeker.init(testPoints);
    QCOMPARE(e, eNoError);

    const QVector<QVector<double>> searchPoint = {{5.0, 5.0}, {5.0, 5.0}};

    const double searchRadius = 2.0;

    QVector<NNSearchResult> result;
    e = seeker.radiusSearch(
            searchPoint,
            searchRadius,
            &result
    );
    QCOMPARE(e, eNoError);

    const int expectedResultSize = 2;
    QCOMPARE(result.size(), expectedResultSize);

    const NNSearchResult &nnSearchResult = result.front();

    std::vector<long> expectedIndexes = {4};
    std::vector<double> expectedValues = {5};
    std::vector<double> expectedDistancesSquared = {0};

    const NNSearchResult &res = result.front();
    QCOMPARE(res.indexes, expectedIndexes);
    QCOMPARE(res.values, expectedValues);
    QCOMPARE(res.distancesSquared, expectedDistancesSquared);

    const double searchRadius2 = 2.01;
    QVector<NNSearchResult> result2;
    e = seeker.radiusSearch(
            searchPoint,
            searchRadius2,
            &result2
    );
    QCOMPARE(e, eNoError);

    const int expectedResultSize2 = 2;
    QCOMPARE(result2.size(), expectedResultSize2);

    std::vector<long> expectedIndexes2 = {4, 3, 5};
    std::vector<double> expectedValues2 = {5, 4, 6};
    std::vector<double> expectedDistancesSquared2 = {0, 2, 2};

    const NNSearchResult &nnSearchResult2 = result2.front();

    QCOMPARE(nnSearchResult2.indexes, expectedIndexes2);
    QCOMPARE(nnSearchResult2.values, expectedValues2);
    QCOMPARE(nnSearchResult2.distancesSquared, expectedDistancesSquared2);
}

void NearestNeighborsSearchTests::readWriteNearestNeighborsTest() {

    ERR_INIT

    const QVector<QPair<double, QVector<double>>> testPoints = buildTestPoints();

    NearestNeighborsSearch seeker;
    e = seeker.init(testPoints);
    QCOMPARE(e, eNoError);

    const QString &testSavePath
            = QDir(qApp->applicationDirPath()).filePath("SoLetItBeWritten.mzml");

    e  = seeker.writeNearestNeighbors(testSavePath);
    QCOMPARE(e, eNoError);

    NearestNeighborsSearch newSeeker;
    e = newSeeker.readNearestNeighbors(
            testSavePath + ".cal",
            testSavePath + ".mat"
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(newSeeker.kdTreeSize(), 10);

}

QTEST_MAIN(NearestNeighborsSearchTests)
#include "NearestNeighborsSearchTests.moc"
