#include "FeatureFinderHillBuilder.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"

#include <QtTest>

#include <iostream>

using namespace Error;


class FeatureFinderHillBuilderTests : public QObject
{
    Q_OBJECT

public:
    FeatureFinderHillBuilderTests() = default;
    ~FeatureFinderHillBuilderTests() override = default;

private slots:
    void buildScanPointGroupsTest();
    void connectCentroidsInGroupedMzValsTest();
    void buildHillsTest();
    void buildHillsRealDataTest();
    void troubleShooting();

private:

    QVector<ScanPoints> testScanPoints();

};


QVector<ScanPoints> FeatureFinderHillBuilderTests::testScanPoints() {

    QVector<ScanPoints> testScanPoints = {
            {ScanPoint(100, 1), ScanPoint(200, 1), ScanPoint(300, 1), ScanPoint(400, 1), ScanPoint(500, 1), ScanPoint(600, 1)},
            {ScanPoint(300, 10), ScanPoint(300.1, 11), ScanPoint(410, 10), ScanPoint(500, 10)},
            {ScanPoint(50, 100), ScanPoint(100, 100), ScanPoint(330, 100), ScanPoint(440, 100)},
            {ScanPoint(44, 1000), ScanPoint(55, 1000), ScanPoint(66, 1000), ScanPoint(100, 1000), ScanPoint(440, 1000), ScanPoint(447, 1000)},
            {ScanPoint(100, 10000), ScanPoint(300, 10000), ScanPoint(447, 10000)}
    };

    return testScanPoints;
}


void FeatureFinderHillBuilderTests::buildScanPointGroupsTest() {

    ERR_INIT

    QVector<ScanPoints> scanPoints = testScanPoints();
    QVector<ScanPoints*> scanPointsPtrs;
    std::transform(
            scanPoints.begin(),
            scanPoints.end(),
            std::back_inserter(scanPointsPtrs),
            [](ScanPoints &sp){return &sp;}
            );

    FeatureFinderParameters params;
    params.skipScanCount = 1;
    params.minScanCount = 2;
    params.tolerancePPM = 1.0;
    params.filterLength = 3;
    params.smoothCount = 1;
    params.sigma = 1;
    params.signalToNoiseRatio = 2;

    FeatureFinderHillBuilder featureFinderHillBuilder;
    e = featureFinderHillBuilder.init(params);
    QCOMPARE(e, eNoError);

    QVector<QVector<QVector<double>>> groupedMzVals;
    QVector<QVector<QVector<float>>> groupedIntensityVals;
    e = featureFinderHillBuilder.buildScanPointGroupsTest(
            scanPointsPtrs,
            &groupedMzVals,
            &groupedIntensityVals
    );
    QCOMPARE(e, eNoError);

    const QVector<int> expectedPointCountPerGroupedScan = {3, 3, 3, 2};

    QCOMPARE(groupedMzVals.size(), expectedPointCountPerGroupedScan.size());
    QCOMPARE(groupedIntensityVals.size(), expectedPointCountPerGroupedScan.size());

    for (int i = 0; i < groupedMzVals.size(); i++) {
        QCOMPARE(groupedMzVals.at(i).size(), expectedPointCountPerGroupedScan.at(i));
    }

    params.skipScanCount = 2;
    e = featureFinderHillBuilder.init(params);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilder.buildScanPointGroupsTest(
            scanPointsPtrs,
            &groupedMzVals,
            &groupedIntensityVals
    );
    QCOMPARE(e, eNoError);

    const QVector<int> expectedPointCountPerGroupedScan2 = {4, 4, 3, 2};

    QCOMPARE(groupedMzVals.size(), expectedPointCountPerGroupedScan2.size());
    QCOMPARE(groupedIntensityVals.size(), expectedPointCountPerGroupedScan2.size());

    for (int i = 0; i < groupedMzVals.size(); i++) {
        QCOMPARE(groupedMzVals.at(i).size(), expectedPointCountPerGroupedScan2.at(i));
    }

    params.skipScanCount = 3;
    e = featureFinderHillBuilder.init(params);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilder.buildScanPointGroupsTest(
            scanPointsPtrs,
            &groupedMzVals,
            &groupedIntensityVals
    );
    QCOMPARE(e, eNoError);

    const QVector<int> expectedPointCountPerGroupedScan3 = {5, 4, 3, 2};

    QCOMPARE(groupedMzVals.size(), expectedPointCountPerGroupedScan3.size());
    QCOMPARE(groupedIntensityVals.size(), expectedPointCountPerGroupedScan3.size());

    for (int i = 0; i < groupedMzVals.size(); i++) {
        QCOMPARE(groupedMzVals.at(i).size(), expectedPointCountPerGroupedScan3.at(i));
    }

    scanPointsPtrs.append(scanPointsPtrs);

    e = featureFinderHillBuilder.buildScanPointGroupsTest(
            scanPointsPtrs,
            &groupedMzVals,
            &groupedIntensityVals
    );
    QCOMPARE(e, eNoError);


}

void FeatureFinderHillBuilderTests::connectCentroidsInGroupedMzValsTest() {

    ERR_INIT

    QVector<ScanPoints> scanPoints = testScanPoints();
    QVector<ScanPoints*> scanPointsPtrs;
    std::transform(
            scanPoints.begin(),
            scanPoints.end(),
            std::back_inserter(scanPointsPtrs),
            [](ScanPoints &sp){return &sp;}
    );

    FeatureFinderParameters params;
    params.skipScanCount = 1;
    params.tolerancePPM = 10.0;
    params.minScanCount = 2;
    params.filterLength = 3;
    params.smoothCount = 1;
    params.sigma = 1;
    params.signalToNoiseRatio = 2;

    FeatureFinderHillBuilder featureFinderHillBuilder;
    e = featureFinderHillBuilder.init(params);
    QCOMPARE(e, eNoError);

    QVector<QVector<QVector<double>>> groupedMzVals;
    QVector<QVector<QVector<float>>> groupedIntensityVals;
    e = featureFinderHillBuilder.buildScanPointGroupsTest(
            scanPointsPtrs,
            &groupedMzVals,
            &groupedIntensityVals
    );
    QCOMPARE(e, eNoError);

    QVector<QVector<int>> connectedCentroidsVecs;
    e = featureFinderHillBuilder.connectCentroidsInGroupedMzValsTest(
            groupedMzVals,
            params.tolerancePPM,
            &connectedCentroidsVecs
    );
    QCOMPARE(e, eNoError);

    const QVector<QVector<int>> expectedResult = {
            {-1, -1, 0, -1, 3, -1, 1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, 3, -1, 4, -1, 0, -1, -1},
            {-1, -1, -1, 0, -1, 2}
    };

    QCOMPARE(connectedCentroidsVecs, expectedResult);
}

void FeatureFinderHillBuilderTests::buildHillsTest() {

    ERR_INIT

    QVector<ScanPoints> scanPoints = testScanPoints();
    QVector<ScanPoints*> scanPointsPtrs;
    std::transform(
            scanPoints.begin(),
            scanPoints.end(),
            std::back_inserter(scanPointsPtrs),
            [](ScanPoints &sp){return &sp;}
    );

    const QList<int> newKeys = {0, 10, 20, 30, 40};

    QMap<int, ScanPoints*> scanNumberVsScanPoints;
    for (int i = 0; i < newKeys.size(); i++) {
        scanNumberVsScanPoints.insert(newKeys.at(i), scanPointsPtrs.at(i));
    }

    FeatureFinderParameters params;
    params.skipScanCount = 2;
    params.tolerancePPM = 100.0;
    params.minScanCount = 2;
    params.filterLength = 3;
    params.smoothCount = 1;
    params.sigma = 1;
    params.signalToNoiseRatio = 2;

    FeatureFinderHillBuilder featureFinderHillBuilder;
    e = featureFinderHillBuilder.init(params);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilder.buildHills(scanNumberVsScanPoints);
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHill*> featureFinderHills;
    e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills);
    QCOMPARE(e, eNoError);
    QCOMPARE(featureFinderHills.size(), 5);

    const auto sortLogic = [](const FeatureFinderHill* l, const FeatureFinderHill *r){
        return l->mzMean() < r->mzMean();
    };

    std::sort(featureFinderHills.begin(), featureFinderHills.end(), sortLogic);

    const QVector<double> expectedMzVals = {100, 300, 440, 447, 500};
    const QVector<QVector<int>> expectedScanNumbersAll = {
            {0, 20, 30, 40},
            {0, 10, 40},
            {20, 30},
            {30, 40},
            {0, 10},
    };

    for (int i = 0; i < featureFinderHills.size(); i++) {

        const FeatureFinderHill* ffh = featureFinderHills.at(i);
        const double expectedMz = expectedMzVals.at(i);
        const QVector<int> &expectedScanNumbers = expectedScanNumbersAll.at(i);

        qDebug() << ffh->mzMean() << ffh->scanNumberIndexes() << ffh->scanNumberIndexMinMax();

        QCOMPARE(ffh->mzMean(), expectedMz);
        QCOMPARE(ffh->scanCount(), expectedScanNumbers.size());
    }

    params.minScanCount = 1;
    FeatureFinderHillBuilder featureFinderHillBuilderSingle;
    e = featureFinderHillBuilderSingle.init(params);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilderSingle.buildHills(scanNumberVsScanPoints);
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHill*> featureFinderHillsSingle;
    e = featureFinderHillBuilderSingle.featureFinderHills(&featureFinderHillsSingle);
    QCOMPARE(e, eNoError);
    QCOMPARE(featureFinderHillsSingle.size(), 15);

    std::sort(featureFinderHillsSingle.begin(), featureFinderHillsSingle.end(), sortLogic);

    for (int i = 0; i < featureFinderHillsSingle.size(); i++) {

        const FeatureFinderHill *ffh = featureFinderHillsSingle.at(i);
        qDebug() << ffh->mzMean() << ffh->scanNumberIndexes() << ffh->scanNumberIndexMinMax();

    }

}

void FeatureFinderHillBuilderTests::buildHillsRealDataTest() {
    QSKIP("TODO: enable with internal test data");

    QSKIP("activate when proper pathing is used");

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    //TODO use proper pathing here
    const QString filePath = "/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq";

    MsReaderParquet msReader;
    msReader.openFile(filePath);

    FeatureFinderParameters params;
    params.skipScanCount = 0;
    params.tolerancePPM = 12.0;
    params.minScanCount = 1;
    params.filterLength = 3;
    params.smoothCount = 1;
    params.sigma = 1;
    params.signalToNoiseRatio = 2;

    FeatureFinderHillBuilder featureFinderHillBuilder;
    e = featureFinderHillBuilder.init(params);
    featureFinderHillBuilder.setRunParallel(false);
    QCOMPARE(e, eNoError);

    const MsLevel msLevel = 1;
    QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints;
    e = msReader.getScanPoints(
            msLevel,
            &scanNumberVsScanPoints
            );
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHill*> featureFinderHills;
    e = featureFinderHillBuilder.buildHills(scanNumberVsScanPoints);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilder.refineHills(true);
    QCOMPARE(e, eNoError);

    e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills);
    QCOMPARE(e, eNoError);

    qDebug() << "Hills found to write" << featureFinderHills.size();

//#define WRITE_TO_MZRT
#ifdef WRITE_TO_MZRT
    e = FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
            msReader.getScanNumberVsScanTime(),
            featureFinderHills,
            "hills.mzrt.csv"
    );
    QCOMPARE(e, eNoError);
#endif

}

void FeatureFinderHillBuilderTests::troubleShooting() {

    ERR_INIT



}


QTEST_MAIN(FeatureFinderHillBuilderTests)

#include "FeatureFinderHillBuilderTests.moc"
