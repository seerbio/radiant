#include "PeakIntegratomatic.h"

#include "ObjectCSVWriters.h"

#include <QHash>
#include <QtTest>

class PeakIntegratomaticTests : public QObject
{
    Q_OBJECT

public:
    PeakIntegratomaticTests() = default;
    ~PeakIntegratomaticTests() override = default;

private slots:

    static void initTest();
    static void simpleIntegratorTest();
    static void simpleIntegratorTest2();

};


void PeakIntegratomaticTests::initTest() {

    ERR_INIT

	PeakIntegratomaticParameters params;
	params.stopThresholdFraction = 1.1;
	params.hysteresis = 0.5;

	PeakIntegratomatic peakIntegratomatic;
	e = peakIntegratomatic.init(params);
	QCOMPARE(e, eValueError);

	params.stopThresholdFraction = 0.1;
	e = peakIntegratomatic.init(params);
	QCOMPARE(e, eNoError);

}

void PeakIntegratomaticTests::simpleIntegratorTest() {

	QSKIP("SDLFKJSD");

    ERR_INIT

	const QVector<int> apexes = {
    	195,206,212,218,226,247,255,264,273,280,294,301,310,320,337,351,357,362,375,379,391,411,425,434,448,457,459,
    	465,504,515,517,540,551,600,607,623,631,664,686,694,751,764,823,827,1145
	};

	const QString signalFilePath = "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/testesProd.csv";
	QVector<float> signal;
	e = ObjectCSVWriters::readVectorFromFile(signalFilePath, &signal);
	QCOMPARE(e, eNoError);

	// const QString signalCountsFilePath = "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/testesCnt.csv";
	// QVector<float> counts;
	// e = ObjectCSVWriters::readVectorFromFile(signalFilePath, &signal);
	// QCOMPARE(e, eNoError);

	PeakIntegratomaticParameters params;
	params.stopThresholdFraction = 0.1;

	PeakIntegratomatic peakIntegratomatic;
	e = peakIntegratomatic.init(params);
	QCOMPARE(e, eNoError);

	QVector<std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex>> peakIntegrationIndexesVsIntensity;
	e = peakIntegratomatic.simpleIntegrator(
		apexes,
		signal.data(),
		signal.size(),
		&peakIntegrationIndexesVsIntensity
		);
	qDebug() << signal.size() << "LSDFJDSL";
	QCOMPARE(e, eNoError);

	for (const std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex> &r : peakIntegrationIndexesVsIntensity) {
		std::cout << "(" << std::get<0>(r).first << "," << std::get<0>(r).second << ")," << std::endl;
	}

}

void PeakIntegratomaticTests::simpleIntegratorTest2() {

	ERR_INIT

	const QString signalFilePath = "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/prod_vec_1078_1211.csv";
	QVector<float> signal;
	e = ObjectCSVWriters::readVectorFromFile(signalFilePath, &signal);
	QCOMPARE(e, eNoError);

    // const QVector<int> apexes = {576,579,608,618,628,646,653,656,663,668,679,694,699,704,709,715,728,742,
    // 	751,757,772,783,787,790,797,801,807,820,826,833,840,846,852,855,858,863,871,887,899,917,928,
    // 	940,947,956,962,968,986,988,997,1004,1008,1034,1067,1104,1115,1138,1156,1166,1178,1200,1222,
    // 	1231,1237,1241,1248,1253,1263,1268,1277,1282,1288,1300,1308,1312,1319,1324,1328,1333,1337,1340
    // 	,1344,1356,1362,1369,1374,1379,1383,1387,1394,1401,1410,1418,1422,1428,1432,1434,1441,1447,1449,1455,
    // 	1459,1462,1467,1474,1480,1486,1493,1498,1502,1509,1513,1519,1525,1532,1540,1544,1552,1558,1561,1564,1581,
    // 	1599,1603,1607,1612,1615,1620,1624,1628,1634,1643,1646,1654,1658,1666,1677,1689,1692,1704,1711,1716,1720,
    // 	1726,1732,1735,1744,1754,1758,1768,1772,1776,1780,1789,1794,1800,1806,1816,1823,1826,1834,1838,1841,1846,
    // 	1850,1854,1858,1868,1876,1882,1887,1893,1914,1935,1957,1975,2000
    // };

	const QVector<int> apexes = {
		1115, 1138, 1156
	};


	PeakIntegratomaticParameters params;
	params.stopThresholdFraction = 0.1;
	params.hysteresis = 0.1;

	PeakIntegratomatic peakIntegratomatic;
	e = peakIntegratomatic.init(params);
	QCOMPARE(e, eNoError);

	QVector<std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex>> peakIntegrationIndexesVsIntensity;
	e = peakIntegratomatic.simpleIntegrator(
		apexes,
		signal.data(),
		signal.size(),
		&peakIntegrationIndexesVsIntensity
		);
	QCOMPARE(e, eNoError);

}


QTEST_MAIN(PeakIntegratomaticTests)

#include "PeakIntegratomaticTests.moc"
