#include "TargetDecoyCandidatePair.h"

#include <QtTest/QtTest>


class TargetDecoyCandidatePairTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairTests() = default;
    ~TargetDecoyCandidatePairTests() override = default;

private Q_SLOTS:

    void gettersTest();
	static void mutateCandidatePeptideTargetTestAccessTest();


};

void TargetDecoyCandidatePairTests::gettersTest() {

    FragLibReaderRow flrr;




	flrr.peptideSequenceChargeKey = "P[Unimod:4]EPTIDE|2";
	flrr.mzVals = {98.06009,227.10268,324.15544,425.20312,538.28718,653.31413,782.35672,49.53370,114.05500,162.58138,213.10522,269.64725,327.16072,391.68202,800.36728,703.31452,574.27193,477.21916,376.17149,263.08742,148.06048,400.68730,352.16092,287.63962,239.11324,188.58940,132.04737,74.53390};
	for (int i = 0; i < flrr.mzVals.size(); i++) {
		flrr.intensityVals.push_back(1.0f/i);
	}

	flrr.ionLabels = {"b1;b2;b3;b4;b5;b6;b7;b1^2;b2^2;b3^2;b4^2;b5^2;b6^2;b7^2;y7;y6;y5;y4;y3;y2;y1;y7^2;y6^2;y5^2;y4^2;y3^2;y2^2;y1^2"};
	flrr.mass = 799.36001	;
	flrr.isDecoy = 0;
	flrr.iRT = 777.0;
	flrr.iM = 666.0;
	flrr.precursorCharge = 2;

	TargetDecoyCandidatePair tdcp;
	tdcp.setFragLibReaderRowPntr(&flrr);
	Err e = tdcp.init();
	QCOMPARE(e, eNoError);
qDebug() << tdcp.peptideStringWithModsDecoy();
    QCOMPARE(tdcp.peptideString(), "PEPTIDE");
	QCOMPARE(tdcp.peptideStringDecoy(), "PDPTIEE");
    QVERIFY(MathUtils::tSame(tdcp.ms2IonsTarget().first().mz, 227.10268f));
    QCOMPARE(tdcp.charge(), 2);
	QVERIFY(MathUtils::tSame(tdcp.mz(false), 400.6873f));
    QVERIFY(MathUtils::tSame(tdcp.mass(false), 799.36001f));
    QVERIFY(MathUtils::tSame(tdcp.iRt(false), 777.0f));
    QCOMPARE(tdcp.totalFragmentCount(), 28);
	QCOMPARE(static_cast<int>(std::abs(tdcp.mass(false) - tdcp.mass(true))), 0);

	QVector<float> ms2IonsExpected = {
		213.087,
		326.171,
		427.219,
		540.303,
		669.345,
		798.388,
		221.121,
		277.663,
		335.176,
		399.698,
		816.399,
		719.346,
		604.319,
		491.235,
		390.187,
		277.103,
		408.703,
		360.177,
		302.663,
		246.121
	};

	const QVector<MS2Ion> ions = tdcp.ms2IonsDecoy();
	QCOMPARE(ions.size(), ms2IonsExpected.size());
	// for (int i = 0; i < ions.size(); i++) {
	// 	QCOMPARE(static_cast<int>(ions[i].mz), static_cast<int>(ms2IonsExpected[i]));
	// }
	QCOMPARE(tdcp.peptideStringDecoy(), "PDPTIEE");
	QCOMPARE(tdcp.peptideStringWithModsDecoy(), "P[Unimod:4]DPTIEE");
	// QVERIFY(MathUtils::tSame(tdcp.mz(true), 408.7029f));
	// QVERIFY(MathUtils::tSame(tdcp.mass(true), 815.39131f));

	//TODO Update the commented out tests when replacing PEPTIDE w/ something that doesn yield an isobaric decoy
}

void TargetDecoyCandidatePairTests::mutateCandidatePeptideTargetTestAccessTest() {

	QVector<MS2Ion> ms2Ions = {
		MS2Ion(845.452, 1, "y7", 1, 1) ,
		MS2Ion(560.304, 0.904366, "y5", 1, 2) ,
		MS2Ion(673.388, 0.676435, "y6", 1, 3) ,
		MS2Ion(423.229, 0.483817, "y7^2", 2, 4) ,
		MS2Ion(473.272, 0.285317, "y4", 1, 5) ,
		MS2Ion(946.499, 0.188432, "y8", 1, 6) ,
		MS2Ion(473.753, 0.155437, "y8^2", 2, 7) ,
		MS2Ion(373.187, 0.130075, "b3", 1, 8) ,
		MS2Ion(486.271, 0.0711044, "b4", 1, 9) ,
		MS2Ion(573.303, 0.0261274, "b5", 1, 10),
		MS2Ion(785.41926, 0.0261274, "b7", 1, 11),
		MS2Ion(899.46219 , 0.0261274, "b8", 1, 12)
	};

	PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("VTWLSPTNK");

	TargetDecoyCandidatePair::mutateCandidatePeptideTargetTestAccess(peptideStringWithMods, ms2Ions);

}


QTEST_MAIN(TargetDecoyCandidatePairTests)
#include "TargetDecoyCandidatePairTests.moc"
