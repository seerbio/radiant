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
	static void trancheMsScanPointsTest();

};


void PrecursorMzArrangetronTests::initTest() {

	QSKIP("Skipping till data is built");

	//TODO build test data

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

	float lastMzPrecursor = -1;
	float lastMzMS2 = -1;
	for (auto it = precursorMzArrangetron.m_precursorMzKeyVsScanNumberMzIntensity.begin();
		it != precursorMzArrangetron.m_precursorMzKeyVsScanNumberMzIntensity.end();
		++it) {

		int key = it.key();
		const QVector<ScanNumberMzIntensity> &v = it.value();

		qDebug()
		<< qPrintable(S_GLOBAL_TIMER.elapsed())
		<< "Precursor Key" << key
		<< "Precursor Mz" << v.front().scanInfoPntr->precursorTargetMz
		<< "Point count" << v.size();

		for (const ScanNumberMzIntensity &snmi : v) {
			const float currentMzPrecursor = snmi.scanInfoPntr->precursorTargetMz;
			const float currentMzMS2 = snmi.mzVal;
			constexpr int precursorMzKeyHashingPrecision = 3;
			e = ErrorUtils::isTrue(
				MathUtils::hashDecimal(lastMzPrecursor, precursorMzKeyHashingPrecision) == MathUtils::hashDecimal(currentMzPrecursor, precursorMzKeyHashingPrecision)
				|| lastMzPrecursor < currentMzPrecursor
				);
			if (e != eNoError) {
				qDebug() << lastMzPrecursor << currentMzPrecursor;
				QCOMPARE(e, eNoError);
			}

			lastMzMS2 = MathUtils::hashDecimal(lastMzPrecursor, precursorMzKeyHashingPrecision)
						!= MathUtils::hashDecimal(currentMzPrecursor, precursorMzKeyHashingPrecision)
					  ? -1
					  : lastMzMS2;

			e = ErrorUtils::isTrue(
			MathUtils::hashDecimal(lastMzMS2, precursorMzKeyHashingPrecision)
					== MathUtils::hashDecimal(currentMzMS2, precursorMzKeyHashingPrecision) || lastMzMS2 < currentMzMS2
				);
			if (e != eNoError) {
				qDebug() << lastMzMS2 << currentMzMS2;
				QCOMPARE(e, eNoError);
			}

			lastMzPrecursor = currentMzPrecursor;
			lastMzMS2 = currentMzMS2;
		}
	}

}

void PrecursorMzArrangetronTests::trancheMsScanPointsTest() {

	// QSKIP("Skipping till data is built");

	//TODO build test data

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

	// constexpr int threadCount = 64;
	// QVector<QVector<ScanNumberMzIntensity*>> scanNumberMzIntensityTranched;
	// e = precursorMzArrangetron.trancheMsScanPoints(threadCount, &scanNumberMzIntensityTranched);
	// QCOMPARE(e, eNoError);

}


QTEST_MAIN(PrecursorMzArrangetronTests)

#include "PrecursorMzArrangetronTests.moc"
