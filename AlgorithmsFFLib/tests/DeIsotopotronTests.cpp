#include <QtTest>
#include <QCoreApplication>

#include "DeIsotopotron.h"
#include "ObjectCSVWriters.h"

class DeIsotopotronTests : public QObject {
Q_OBJECT

public:
    DeIsotopotronTests() = default;

    ~DeIsotopotronTests() override = default;

private slots:

    static void deisotopeTest();

};

void DeIsotopotronTests::deisotopeTest() {

    ERR_INIT

	const QString scanPointsFile
		= "/home/andrewnichols/Repos/Builds/PythiaDIA/bin/Scan27355.txt"; //TODO put this in tests

	ScanPoints scanPoints;
	e = ObjectCSVWriters::readScanPoints(scanPointsFile, &scanPoints);
	QCOMPARE(e, eNoError);
	QCOMPARE(scanPoints.size(), 91);

	constexpr float ppmTol = 10.0f;
	e = DeIsotopotron::deisotopeTandemScan(ppmTol, &scanPoints);
	QCOMPARE(e, eNoError);


}


QTEST_MAIN(DeIsotopotronTests)

#include "DeIsotopotronTests.moc"
