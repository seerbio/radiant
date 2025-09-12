//
// Created by anichols on 11/07/2021.
//

#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertronV2.h"
#include "MsCalibratomatic.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"


class MsCalibratomaticSettertronV2Tests : public QObject
{
    Q_OBJECT

public:
    MsCalibratomaticSettertronV2Tests() = default;
    ~MsCalibratomaticSettertronV2Tests() override = default;

private Q_SLOTS:

void testme();

};

void MsCalibratomaticSettertronV2Tests::testme() {

	ERR_INIT

	const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
	// const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";

	MsReaderPointerAcc reader;
	e = reader.openFile(filename);
	QCOMPARE(e, eNoError);

	// const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/human_plasma_arath_entrapment.fasta.predicted.speclib";
	const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF";
	QVector<FragLibReaderRow> fragLibReaderRows;
	e = FragLibReader::getFragLibReaderRows(
		fragLibUri,
		&fragLibReaderRows
		);
	QCOMPARE(e, eNoError);

	PythiaParameters pythiaParameters = PythiaParameterReader::genericPythiaParametersForTests();
	// pythiaParameters.ms2ExtractionWidthPPM = 8;

	TargetDecoyCandidatePairManager tdcpManager;
	e = tdcpManager.init(
		pythiaParameters,
		&fragLibReaderRows
		);
	QCOMPARE(e, eNoError);

	constexpr int msLevel = 1;
	const QMap<ScanNumber, MsScanInfo*> ms1ScanInfos = reader.ptr->getMsScanInfos(msLevel);
	QVector<MsScan> msScans;

	e = reader.ptr->extractScanPoints(
		ms1ScanInfos.values().toVector(),
		&msScans
		);
	QCOMPARE(e, eNoError);

	MsFrameV2 msFrameMS1;
	e = msFrameMS1.init(msScans);
	QCOMPARE(e, eNoError);


	// MsCalibratomaticSettertronV2 setter;
	// e = setter.init(
	// 	&tdcpManager,
	// 	&reader,
	// 	&pythiaParameters,
	// 	&msFrameMS1
	// 	);
	// QCOMPARE(e, eNoError);
	//
	// MsCalibratomatic msCalibratomatic;
	// e = setter.buildMsCalibratomatic(&msCalibratomatic);
	// QCOMPARE(e, eNoError);


}

QTEST_MAIN(MsCalibratomaticSettertronV2Tests)
#include "MsCalibratomaticSettertronV2Tests.moc"
