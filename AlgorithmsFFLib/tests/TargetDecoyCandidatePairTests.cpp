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
	static void mutateCandidatePeptideTargetTestAccessTestPenultimate();
	static void mutateCandidatePeptideTargetTestAccessTestTerminalByPenultimate();


};

void TargetDecoyCandidatePairTests::gettersTest() {

    QSKIP("fix test");


    MS2Ion ms2Ion1;
    ms2Ion1.mz = 666.6;

    MS2Ion ms2Ion2;
    ms2Ion2.mz = 777.7;

    const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods("CHAUNCYANDFLOPS");
    const QVector<MS2Ion> ms2IonsTarget = {ms2Ion1};
    const QVector<MS2Ion> ms2IonsDecoy = {ms2Ion2};
    const int charge = 2;
    const double mass = 666.6;
    const double iRt = 66.6;
    const int totalFramentCount = 12;

    // TargetDecoyCandidatePair targetDecoyCandidatePair(
    //         peptideStringWithMods,
    //         ms2IonsTarget,
    //         ms2IonsDecoy,
    //         charge,
    //         mass,
    //         iRt,
    //         iRt,
    //         totalFramentCount,
    //         0.0
    //         );
    //
    // QCOMPARE(targetDecoyCandidatePair.peptideStringWithMods(), peptideStringWithMods);
    // QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsTarget().first().mz, 666.6f));
    // QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.ms2IonsDecoy().first().mz, 777.7f));
    // QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mz(false), 334.3072f));
    // QCOMPARE(targetDecoyCandidatePair.charge(), charge);
    // QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.mass(), 666.6f));
    // QVERIFY(MathUtils::tSame(targetDecoyCandidatePair.iRt(), 66.6f));
    // QCOMPARE(targetDecoyCandidatePair.totalFragmentCount(), 12);

}

void TargetDecoyCandidatePairTests::mutateCandidatePeptideTargetTestAccessTestPenultimate() {
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

	DecoyFragmentShiftMode shiftMode = DecoyFragmentShiftMode::ShiftPenultimate;

    const QVector<MS2Ion> decoys = TargetDecoyCandidatePair::mutateCandidatePeptideTargetTestAccess(
            peptideStringWithMods,
            ms2Ions,
            shiftMode
            );
    QCOMPARE(decoys.size(), ms2Ions.size());

    constexpr float tolerance = 0.001f;
    const auto assertMzNear = [tolerance](const float observed, const float expected) {
        QVERIFY2(
                std::abs(observed - expected) <= tolerance,
                qPrintable(QString("Observed m/z %1 differs from expected %2").arg(observed).arg(expected))
                );
    };

    const QVector<float> expectedDecoyMzs = {
            859.46765f, // y7
            574.31965f, // y5
            687.40365f, // y6
            430.23682f, // y7^2
            487.28765f, // y4
            946.49900f, // y8
            473.75300f, // y8^2
            359.17135f, // b3
            472.25535f, // b4
            559.28735f, // b5
            771.40361f, // b7
            899.46219f  // b8
    };
    QCOMPARE(expectedDecoyMzs.size(), decoys.size());
    for (int i = 0; i < decoys.size(); ++i) {
        assertMzNear(decoys.at(i).mz, expectedDecoyMzs.at(i));
    }

}

void TargetDecoyCandidatePairTests::mutateCandidatePeptideTargetTestAccessTestTerminalByPenultimate() {
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

	DecoyFragmentShiftMode shiftMode = DecoyFragmentShiftMode::ShiftTerminalByPenultimate;

	const QVector<MS2Ion> decoys = TargetDecoyCandidatePair::mutateCandidatePeptideTargetTestAccess(
			peptideStringWithMods,
			ms2Ions,
			shiftMode
			);
	QCOMPARE(decoys.size(), ms2Ions.size());

	constexpr float tolerance = 0.001f;
	const auto assertMzNear = [tolerance](const float observed, const float expected) {
		QVERIFY2(
				std::abs(observed - expected) <= tolerance,
				qPrintable(QString("Observed m/z %1 differs from expected %2").arg(observed).arg(expected))
				);
	};

	const QVector<float> expectedDecoyMzs = {
		845.45200f, // y7
		574.31965f, // y5
		687.40365f, // y6
		423.22900f, // y7^2
		487.28765f, // y4
		946.49900f, // y8
		473.75300f, // y8^2
		359.17135f, // b3
		472.25535f, // b4
		559.28735f, // b5
		785.41926f, // b7
		899.46219f  // b8
};
	QCOMPARE(expectedDecoyMzs.size(), decoys.size());
	for (int i = 0; i < decoys.size(); ++i) {
		assertMzNear(decoys.at(i).mz, expectedDecoyMzs.at(i));
	}

}


QTEST_MAIN(TargetDecoyCandidatePairTests)
#include "TargetDecoyCandidatePairTests.moc"
