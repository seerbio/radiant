//
// Created by anichols on 11/07/2021.
//

#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertronV2.h"
#include "MsReaderPointerAcc.h"


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

	MsReaderPointerAcc reader;
	e = reader.openFile(filename);
	QCOMPARE(e, eNoError);

	const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/human_plasma_arath_entrapment.fasta.predicted.speclib";
	QVector<FragLibReaderRow> fragLibReaderRows;
	e = FragLibReader::getFragLibReaderRows(
		fragLibUri,
		&fragLibReaderRows
		);
	QCOMPARE(e, eNoError);

	MsCalibratomaticSettertronV2 setter;
	e = setter.init(&reader);
	QCOMPARE(e, eNoError);


}

QTEST_MAIN(MsCalibratomaticSettertronV2Tests)
#include "MsCalibratomaticSettertronV2Tests.moc"
