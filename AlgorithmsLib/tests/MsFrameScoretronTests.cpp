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

private:

    static PythiaParameters pythiaParameters() {

        const QString &paramsFile
                = QDir(qApp->applicationDirPath()).filePath("WorkFlowTestsParams.pythia");

        PythiaParameterReader reader;
        PythiaParameters pythiaParameters;
        reader.readFile(paramsFile);
        reader.loadPythiaParameters(&pythiaParameters);

        pythiaParameters.topNMs2Ions = 12;
        pythiaParameters.print();
        return pythiaParameters;
    }

};

void MsFrameScoretronTests::scoreCandidatesTest() {

    ERR_INIT

    const QString mzMLFileURI
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/RawData/2022_02_22_Homo_sapiens_UP000005640.fasta.fragLib");

    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "474966";

    MsFrameScoretron msFrameScoretron;
    QPair<Err, QVector<PSMsReaderRow>> result = msFrameScoretron.scoreCandidates(
            pythiaParameters(),
            mzMLFileURI,
            fragLibPath,
            uniqueMsInfoScanKey,
            {474.966 - 5.5, 474.966 + 5.5},
            false
            );
    QCOMPARE(result.first, eNoError);

}


QTEST_MAIN(MsFrameScoretronTests)

#include "MsFrameScoretronTests.moc"
