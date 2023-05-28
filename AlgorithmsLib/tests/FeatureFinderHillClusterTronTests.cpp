#include "FeatureFinderHillClusterTron.h"
#include "FeatureFinderHillBuilder.h"


#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"

#include <QtTest>

#include <iostream>

using namespace Error;


class FeatureFinderHillClusterTronTests : public QObject
{
    Q_OBJECT

public:
    FeatureFinderHillClusterTronTests() = default;
    ~FeatureFinderHillClusterTronTests() override = default;

private slots:

    void clusterHillsToyDataTest();
    void clusterHillsRealDataTest();


};

namespace {

    FeatureFinderHill buildHill(
            double mz,
            double intensity,
            QPair<int, int> startStopFrameIndex
            ) {

        const int hillLength = startStopFrameIndex.second -startStopFrameIndex.first;
        const int hillLengthMidPoint = std::floor(hillLength / 2.0);

        const double intensityIncrement = intensity / hillLengthMidPoint;

        int counter = 0;

        FeatureFinderHill featureFinderHill;

        for (int i = startStopFrameIndex.first; i <= startStopFrameIndex.second; i++) {

            const int corrector = hillLength % 2 == 0 ? 2 : 3;

            const double intensityInc = counter++ <= hillLengthMidPoint
                    ? counter * intensityIncrement
                    : intensity - ((counter - hillLengthMidPoint - corrector) * intensityIncrement);

            featureFinderHill.addPoint(
                    i,
                    i,
                    mz,
                    intensityInc
                    );
        }


        return featureFinderHill;
    }

}//namespace
void FeatureFinderHillClusterTronTests::clusterHillsToyDataTest() {

    const FeatureFinderHill featureFinderHill1 = buildHill(
            100.01,
            1000.0,
            {10, 20}
            );

    const FeatureFinderHill featureFinderHill2 = buildHill(
            101.01,
            500.0,
            {12, 18}
    );

    const FeatureFinderHill featureFinderHill3 = buildHill(
            102.01,
            250.0,
            {13, 17}
    );

    const FeatureFinderHill featureFinderHillInerferer = buildHill(
            100.51,
            500.0,
            {14, 22}
    );

    const QVector<FeatureFinderHill> featureFinderHills = {
            featureFinderHill1,
            featureFinderHill2,
            featureFinderHill3,
            featureFinderHillInerferer
    };

    ERR_INIT

    FeatureFinderParameters featureFinderParameters;
    featureFinderParameters.tolerancePPM = 1000.0;
    featureFinderParameters.skipScanCount = 0;
    featureFinderParameters.minScanCount = 1;
    featureFinderParameters.filterLength = 2;
    featureFinderParameters.smoothCount = 1;
    featureFinderParameters.sigma = 1;
    featureFinderParameters.signalToNoiseRatio = 1.0;
    featureFinderParameters.cosineSimThreshold = 0.8;

    FeatureFinderHillClusterTron featureFinderHillClusterTron;
    e = featureFinderHillClusterTron.init(featureFinderParameters);
    QCOMPARE(e, eNoError);

    QVector<HillsClustering> hillClustersByIndexs;
    QMap<int, FeatureFinderHill> featureFinderHillsMap;
    e = featureFinderHillClusterTron.clusterHills(
            featureFinderHills,
            true,
            &hillClustersByIndexs,
            &featureFinderHillsMap
            );
    QCOMPARE(e, eNoError);

    QCOMPARE(hillClustersByIndexs.size(), 2);

    const HillsClustering &firstRes = hillClustersByIndexs.front();
    qDebug() << firstRes.hillsMapIndexes;
    QCOMPARE(firstRes.hillsMapIndexes.size(), 3);
    QCOMPARE(firstRes.hillsMapIndexes.front(), 0);
    QCOMPARE(firstRes.hillsMapIndexes.at(1), 1);
    QCOMPARE(firstRes.hillsMapIndexes.back(), 2);
    QCOMPARE(firstRes.cosineSim,  0.8669);
    QCOMPARE(firstRes.monoOffset, 0);
    QCOMPARE(firstRes.charge, 1);

    const HillsClustering &secondRes = hillClustersByIndexs.back();
    QCOMPARE(secondRes.hillsMapIndexes.size(), 1);
    QCOMPARE(secondRes.hillsMapIndexes.front(), 3);
    QCOMPARE(secondRes.cosineSim, -1.0);
    QCOMPARE(secondRes.monoOffset, -1);
    QCOMPARE(secondRes.charge, -1);

}


void FeatureFinderHillClusterTronTests::clusterHillsRealDataTest() {

    ERR_INIT

    QSKIP("later");

    //TODO use proper pathing here
    const QString filePath = "/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq";

    MsReaderParquet msReader;
    msReader.openFile(filePath);

    FeatureFinderParameters params;
    params.skipScanCount = 0;
    params.tolerancePPM = 12.0;
    params.minScanCount = 1;
//    params.signalToNoiseRatio = 1;

    FeatureFinderHillBuilder featureFinderHillBuilder;
    e = featureFinderHillBuilder.init(params);
    featureFinderHillBuilder.setRunParallel(false);
    QCOMPARE(e, eNoError);

    const MsLevel msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    e = msReader.getScanPoints(
            msLevel,
            &scanNumberVsScanPoints
            );
    QCOMPARE(e, eNoError);

    QVector<FeatureFinderHill> featureFinderHills;
    e = featureFinderHillBuilder.buildHills(
            scanNumberVsScanPoints,
            &featureFinderHills
            );
    QCOMPARE(e, eNoError);

    qDebug() << "Hills found to write" << featureFinderHills.size();

#define WRITE_TO_MZRT
#ifdef WRITE_TO_MZRT
    e = FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
            msReader.getScanNumberVsScanTime(),
            featureFinderHills,
            "hills.mzrt.csv"
    );
    QCOMPARE(e, eNoError);
#endif

}


QTEST_MAIN(FeatureFinderHillClusterTronTests)

#include "FeatureFinderHillClusterTronTests.moc"
