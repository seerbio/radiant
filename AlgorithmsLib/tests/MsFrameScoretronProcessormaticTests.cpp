#include "MsFrameScoretronProcessormatic.h"

#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoreVectorReader.h"
#include "MsReaderParquet.h"
#include "ParquetReader.h"

#include <QElapsedTimer>
#include <QVector>
#include <QtTest>


class MsFrameScoretronProcessormaticTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronProcessormaticTests() = default;
    ~MsFrameScoretronProcessormaticTests() override = default;


private slots:

    void rescoreMsFrameTest();


};


void MsFrameScoretronProcessormaticTests::rescoreMsFrameTest() {

    ERR_INIT

    const QString scoreVectorsFilePath
            = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.reCal.805116.frameScores";

    const QString extractsFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq.reCal.805116.extracts";

    const QString msDataFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq";

    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "805116";
    const QPair<double, double> mzTargetStartStop = {804.0, 806.0};

    MsFrameScoretronProcessormatic msFrameScoretronProcessormatic;
    msFrameScoretronProcessormatic.init(
            scoreVectorsFilePath,
            extractsFilePath,
            PythiaParameterReader::genericPythiaParametersForTests(),
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop
            );
    QCOMPARE(e, eNoError);

    QVector<PSMsReaderRow> psmReaderRows;
    e = msFrameScoretronProcessormatic.processFrameScoreVectors(&psmReaderRows);
    QCOMPARE(e, eNoError);



}


QTEST_MAIN(MsFrameScoretronProcessormaticTests)

#include "MsFrameScoretronProcessormaticTests.moc"
