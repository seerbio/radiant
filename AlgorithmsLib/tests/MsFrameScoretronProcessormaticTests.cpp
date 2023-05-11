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

    const QString fragLibFilePath = "/home/anichols/Desktop/2022_02_22_Homo_sapiens_UP000005640.fragLib";
    const QString msDataFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq";
    const double scanTime = 20.5298;

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    const ScanNumber scanNumber = msReaderParquet.getNearestScanNumberFromScanTime(scanTime);

    MsScanInfo msScanInfo;
    e = msReaderParquet.getMsScanInfo(scanNumber, &msScanInfo);
    QCOMPARE(e, eNoError);

    qDebug() << "TargetKey" << msScanInfo.targetScanKey();

    const PythiaParameters params = PythiaParameterReader::genericPythiaParametersForTests();
    e = ErrorUtils::isTrue(params.isValid());
    QCOMPARE(e, eNoError);
    params.print();

    const QPair<double, double> mzTargetStartStop = {
            msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower,
            msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper
    };

    QMap<PeptideStringWithMods, bool> fragPredsIsDecoy;
    QMap<PeptideStringWithMods, QVector<MS2Ion>> fragPreds;
    e = MsFrameScoretron::buildFragIonLibForTargetMz(
            params,
            fragLibFilePath,
            mzTargetStartStop,
            &fragPreds,
            &fragPredsIsDecoy
    );
    QCOMPARE(e, eNoError);

    const bool applySmooth2D = true;
    const QString uniqueMsInfoScanKey = msScanInfo.targetScanKey();

    MsFrame msFrame;
    e = MsFrameScoretron::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            params,
            mzTargetStartStop,
            applySmooth2D,
            &msFrame
    );
    QCOMPARE(e, eNoError);

    MsFrameScoretronProcessormatic msFrameScoretronProcessormatic;
            msFrameScoretronProcessormatic.init(
                    fragPreds,
                    msFrame,
                    params,
                    scoreVectorsFilePath
                    );
    QCOMPARE(e, eNoError);

    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> topCansInFrameIndex;
    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>>> topCansInFrameIndexVsDiscScore;
    msFrameScoretronProcessormatic.processLogicForFrameScores(
            &topCansInFrameIndex,
            &topCansInFrameIndexVsDiscScore
            );
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(MsFrameScoretronProcessormaticTests)

#include "MsFrameScoretronProcessormaticTests.moc"
