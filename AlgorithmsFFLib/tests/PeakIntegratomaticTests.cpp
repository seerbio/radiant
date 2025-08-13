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

	QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexesVsIntensity;
	e = peakIntegratomatic.simpleIntegrator(
		apexes,
		signal.data(),
		signal.size(),
		&peakIntegrationIndexesVsIntensity
		);
	qDebug() << signal.size() << "LSDFJDSL";
	QCOMPARE(e, eNoError);

	for (const QPair<PeakIntegrationIndexes, float> &r : peakIntegrationIndexesVsIntensity) {
		std::cout << "(" << r.first.first << "," << r.first.second << ")," << std::endl;
	}



}


QTEST_MAIN(PeakIntegratomaticTests)

#include "PeakIntegratomaticTests.moc"
