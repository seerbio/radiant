#include "Error.h"
#include "MS2PointsExtractomatic.h"

#include <QtTest>


class MS2PointsExtractomaticTests : public QObject
{
    Q_OBJECT

public:
    MS2PointsExtractomaticTests() = default;
    ~MS2PointsExtractomaticTests() override = default;

private slots:

    void initTest();
    void findScanPointsTest();
    void deletePointsTest();
    void updateScanPointsTest();

private:

    const QVector<ScanPoints> m_scanPoints = {
            {{1000.1, 1000}, {1001.1, 2000}},
            {{100.1, 10000}, {101.1, 20000}}
    };

};


void MS2PointsExtractomaticTests::initTest() {

    ERR_INIT

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init(m_scanPoints);
    QCOMPARE(e, eNoError);

    e = ms2PointsExtractomatic.init({});
    QCOMPARE(e, eEmptyContainerError);

}

void MS2PointsExtractomaticTests::findScanPointsTest() {

    ERR_INIT

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init(m_scanPoints);
    QCOMPARE(e, eNoError);

    QMap<Id, ScanPoint> foundIdVsScanPoints;
    e = ms2PointsExtractomatic.findScanPoints(
            1000.0,
            1000.101,
            &foundIdVsScanPoints
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(foundIdVsScanPoints.firstKey(), 0);
    QCOMPARE(foundIdVsScanPoints.first().y(), 1000.0);

}

void MS2PointsExtractomaticTests::deletePointsTest() {

    ERR_INIT

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init(m_scanPoints);
    QCOMPARE(e, eNoError);

    QMap<Id, ScanPoint> foundIdVsScanPoints;
    e = ms2PointsExtractomatic.findScanPoints(
            1000.0,
            1000.101,
            &foundIdVsScanPoints
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(foundIdVsScanPoints.firstKey(), 0);
    QCOMPARE(foundIdVsScanPoints.first().y(), 1000.0);

    e = ms2PointsExtractomatic.removeScanPoints(foundIdVsScanPoints);
    QCOMPARE(e, eNoError);

}

void MS2PointsExtractomaticTests::updateScanPointsTest() {

    ERR_INIT

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init(m_scanPoints);
    QCOMPARE(e, eNoError);

    QMap<Id, ScanPoint> foundIdVsScanPoints;
    e = ms2PointsExtractomatic.findScanPoints(
            1000.0,
            1000.101,
            &foundIdVsScanPoints
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(foundIdVsScanPoints.firstKey(), 0);
    QCOMPARE(foundIdVsScanPoints.first().y(), 1000.0);

    foundIdVsScanPoints[0].ry() = 666.0;

    e = ms2PointsExtractomatic.updateScanPoints(foundIdVsScanPoints);
    QCOMPARE(e, eNoError);

    foundIdVsScanPoints.clear();
    e = ms2PointsExtractomatic.findScanPoints(
            1000.0,
            1000.101,
            &foundIdVsScanPoints
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(foundIdVsScanPoints.firstKey(), 0);
    QCOMPARE(foundIdVsScanPoints.first().y(), 666.0);

}


QTEST_MAIN(MS2PointsExtractomaticTests)

#include "MS2PointsExtractomaticTests.moc"
