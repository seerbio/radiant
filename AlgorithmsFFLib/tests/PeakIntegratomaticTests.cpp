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

	const QString signalFilePath = "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/testes.csv";
	QVector<float> signal;
	e = ObjectCSVWriters::readVectorFromFile(signalFilePath, &signal);
	QCOMPARE(e, eNoError);

	PeakIntegratomaticParameters params;

	PeakIntegratomatic peakIntegratomatic;
	e = peakIntegratomatic.init(params);
	QCOMPARE(e, eNoError);

	QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexesVsIntensity;
	e = peakIntegratomatic.simpleIntegrator(signal.data(), signal.size(), &peakIntegrationIndexesVsIntensity);
	qDebug() << signal.size() << "LSDFJDSL";
	QCOMPARE(e, eNoError);

	for (const QPair<PeakIntegrationIndexes, float> &r : peakIntegrationIndexesVsIntensity) {
		std::cout << "(" << r.first.first << "," << r.first.second << ")," << std::endl;
	}



}


QTEST_MAIN(PeakIntegratomaticTests)

#include "PeakIntegratomaticTests.moc"
