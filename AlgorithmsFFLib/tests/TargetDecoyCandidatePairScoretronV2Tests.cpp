#include "TargetDecoyCandidatePairScoretronV2.h"

#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"

#include <QtTest/QtTest>

#include <iostream>

class TargetDecoyCandidatePairScoretronV2Tests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairScoretronV2Tests() = default;
    ~TargetDecoyCandidatePairScoretronV2Tests() override = default;

private Q_SLOTS:
    static void buildIntegrationVecCosineSimTest();


};

void TargetDecoyCandidatePairScoretronV2Tests::buildIntegrationVecCosineSimTest() {

	ERR_INIT

	// QSKIP("Replace file w/ smaller one");

	const QString filename
		= "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";

	MsReaderPointerAcc reader;
	e = reader.openFile(filename);
	QCOMPARE(e, eNoError);


	QMap<ScanNumber, MsScanInfo> msScanInfos = reader.ptr->getMsScanInfos(2);

	QMap<MzTargetKey, QVector<MsScanInfo*>> mzTargetKeyVsMsScanInfos;
	for (MsScanInfo &msi : msScanInfos) {
		mzTargetKeyVsMsScanInfos[msi.targetKey()].push_back(&msi);
	}

	TargetDecoyCandidatePairScoretronV2 tdcpScoretron;
	e = tdcpScoretron.init(
		mzTargetKeyVsMsScanInfos,
		PythiaParameterReader::genericPythiaParametersForTests(),
		16,
		&reader
		);
	QCOMPARE(e, eNoError);

	for (int i = 0; i < 16; i ++) {
		QVector<float> v(24, static_cast<float>(i));
		std::copy_n(v.begin(), 24, tdcpScoretron.m_xicsAlignasIntensity[i]);
	}

	tdcpScoretron.m_targetFrameIndexMax = 24;
	tdcpScoretron.m_ms2IonsCount = 16;

	// e = tdcpScoretron.buildIntegrationVecCosineSim();
	QCOMPARE(e, eNoError);

}


QTEST_MAIN(TargetDecoyCandidatePairScoretronV2Tests)
#include "TargetDecoyCandidatePairScoretronV2Tests.moc"
