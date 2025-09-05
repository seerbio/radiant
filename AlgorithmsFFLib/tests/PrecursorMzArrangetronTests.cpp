#include <QtTest>
#include <QCoreApplication>

#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
#include "PrecursorMzArrangetron.h"

class PrecursorMzArrangetronTests : public QObject {
Q_OBJECT

public:
    PrecursorMzArrangetronTests() = default;

    ~PrecursorMzArrangetronTests() override = default;

private slots:

    static void initTest();

};

void PrecursorMzArrangetronTests::initTest() {

    ERR_INIT

	QString msFile = "/home/andrewnichols/Desktop/Data/MsData/EXP21126_2021ms0425XRC6_C.raw.mzML";
	msFile = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
	// msFile = "/home/andrewnichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML";


	MsReaderPointerAcc msReader;
	e = msReader.openFile(msFile);
	QCOMPARE(e, eNoError);

	PrecursorMzArrangetron precursorMzArrangetron;
	e = precursorMzArrangetron.init(&msReader);
	QCOMPARE(e, eNoError);

}


QTEST_MAIN(PrecursorMzArrangetronTests)

#include "PrecursorMzArrangetronTests.moc"
