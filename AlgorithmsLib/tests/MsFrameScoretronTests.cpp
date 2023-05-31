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

    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Desktop/Testing/2022_02_22_Homo_sapiens_UP000005640.fragLib");
            = QStringLiteral("/home/anichols/Downloads/2022-04-23-decoys-contam-UP000005640_9606.target.20210602.human_plasma.fasta.csv.025.fragLib");

    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "635039";

    MsFrameScoretron msFrameScoretron;
    e = msFrameScoretron.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            msDataFilePath,
            fragLibPath,
            uniqueMsInfoScanKey,
            {635.039 - 5.5, 635.039 + 5.5}
            );
    QCOMPARE(e, eNoError);

    QString frameScoreVectorsFilePath;
    e = msFrameScoretron.extractHillsForCandidtates(&frameScoreVectorsFilePath);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(MsFrameScoretronTests)

#include "MsFrameScoretronTests.moc"
