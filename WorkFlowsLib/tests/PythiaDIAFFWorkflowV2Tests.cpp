#include "PythiaDIAFFWorkflowV2.h"

#include "Error.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"
#include "PythiaDIAFFWorkflowV2.h"

#include <QtTest>

class PythiaDIAFFWorkflowV2Tests : public QObject
{
    Q_OBJECT

public:

    PythiaDIAFFWorkflowV2Tests() = default;
    ~PythiaDIAFFWorkflowV2Tests() override = default;

private slots:

    static void initTest();

    void buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest();


};

void PythiaDIAFFWorkflowV2Tests::initTest() {
    ERR_INIT

	const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
	// const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";

	const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF";
	// const QString fragLibUri  = "/home/andrewnichols/Desktop/Data/Libraries/human_plasma_arath_entrapment.fasta.predicted.speclib";

	const PythiaParameters pythiaParameters = PythiaParameterReader::genericPythiaParametersForTests();

	PythiaDIAFFWorkflowV2 pythiaDIAFFWorkflowV2;
	e = pythiaDIAFFWorkflowV2.init(pythiaParameters, fragLibUri);
	QCOMPARE(e, eNoError);

	e = pythiaDIAFFWorkflowV2.processFile(filename);
	QCOMPARE(e, eNoError);


}

void PythiaDIAFFWorkflowV2Tests::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest() {}


QTEST_MAIN(PythiaDIAFFWorkflowV2Tests)

#include "PythiaDIAFFWorkflowV2Tests.moc"
