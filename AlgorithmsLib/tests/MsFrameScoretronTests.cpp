#include "MsFrameScoretron.h"

#include <QElapsedTimer>
#include <QtTest>


class MsFrameScoretronTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronTests() = default;
    ~MsFrameScoretronTests() override = default;


private slots:

    void scoreCandidatesTest();

};

void MsFrameScoretronTests::scoreCandidatesTest() {

    ERR_INIT

    const QString msDataFilePath
            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");

    const QString fragLibBackgroundPath
            = QStringLiteral("/home/anichols/Desktop/Libraries/uniparc_upid_UP000027126_2023_07_02.fasta.csv.fragLib");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");

//    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "444952";
//    const double target = 444.952;

    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "454957";
    const double target = 454.957;

    const double targetWindowSize = 5.5;

    MsFrameScoretron msFrameScoretron;
    e = msFrameScoretron.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            msDataFilePath,
            fragLibPath,
            fragLibBackgroundPath,
            uniqueMsInfoScanKey,
            {target - targetWindowSize, target + targetWindowSize}
            );
    QCOMPARE(e, eNoError);

    QMap<FrameIndex , QVector<ScoredCandidate>> frameIndexVsScoredCandidates;
    e = msFrameScoretron.scoreFrameCandidates(&frameIndexVsScoredCandidates);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(MsFrameScoretronTests)

#include "MsFrameScoretronTests.moc"
