#include "PythiaDDAWorkflow.h"

#include "Error.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

#include <QtTest>

class PythiaDDAWorkflowTests : public QObject
{
    Q_OBJECT

public:

    PythiaDDAWorkflowTests() = default;
    ~PythiaDDAWorkflowTests() override = default;

private slots:

    static void initTest();
    static void processFileTest();



};

void PythiaDDAWorkflowTests::initTest() {

	ERR_INIT




}

void PythiaDDAWorkflowTests::processFileTest() {

	ERR_INIT

	// const QString ddaMsDataFilePath = "/home/andrewnichols/Desktop/Data/MsData/EXP21126_2021ms0425XRC6_C.raw.mzML";
	const QString ddaMsDataFilePath = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
	const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF";

	const PythiaParameters parameters = PythiaParameterReader::genericPythiaParametersForTests();

	PythiaDDAWorkflow ddaMsWorkflow;
	e = ddaMsWorkflow.init(
		parameters,
		fragLibUri
		);
	QCOMPARE(e, eNoError);

	e = ddaMsWorkflow.processFile(ddaMsDataFilePath);
	QCOMPARE(e, eNoError);

}


QTEST_MAIN(PythiaDDAWorkflowTests)

#include "PythiaDDAWorkflowTests.moc"
