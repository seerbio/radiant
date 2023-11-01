#include "Error.h"
#include "GlobalSettings.h"
#include "ConvertMzMLToParquetWorkFlow.h"

#include <QtTest>


using namespace Error;


class ConvertMzMLToParquetWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    ConvertMzMLToParquetWorkFlowTests() = default;
    ~ConvertMzMLToParquetWorkFlowTests() override = default;

private slots:

    void convertMzMLToParquetRunTest();

};

void ConvertMzMLToParquetWorkFlowTests::convertMzMLToParquetRunTest() {

    ERR_INIT

    QSKIP("activate when proper pathing is used");

    //TODO use proper path procedures after finding small file.
    const QString mzMLFilepath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML");

    const QString expectedOutputFilePath
            = mzMLFilepath + S_GLOBAL_SETTINGS.DOT_PRQ;

//    const QString &fastaFilePath
//            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

//    QString outputFilePath;
//    e = ConvertMzMLToParquetWorkFlow::convertMzMLToParquet(
//            mzMLFilepath,
//            &outputFilePath
//            );
//    QCOMPARE(e, eNoError);
//    QCOMPARE(outputFilePath, expectedOutputFilePath);
//
//    e = ErrorUtils::fileExists(expectedOutputFilePath);
//    QCOMPARE(e, eNoError);

}


QTEST_MAIN(ConvertMzMLToParquetWorkFlowTests)

#include "ConvertMzMLToParquetWorkFlowTests.moc"
