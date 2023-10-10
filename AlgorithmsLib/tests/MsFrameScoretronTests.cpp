#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "FragLibReader.h"
#include "MsReaderParquet.h"

#include <QElapsedTimer>
#include <QtTest>

#include <iostream>

class MsFrameScoretronTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronTests() = default;
    ~MsFrameScoretronTests() override = default;


private slots:

    void scoreCandidatesTest();
    void scoreCandidatesRecalTest();

};

void MsFrameScoretronTests::scoreCandidatesTest() {

    QSKIP("undo me baby");

    ERR_INIT

//    const QString msDataFilePath
//            = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");
//
//    const QString fragLibBackgroundPath
//            = QStringLiteral("/home/anichols/Desktop/Libraries/uniparc_upid_UP000027126_2023_07_02.fasta.csv.fragLib");
//
//    const QString fragLibPath
//            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");
//
////    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "444952";
////    const double target = 444.952;
//
//    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "454957";
//    const double target = 454.957;
//
//    const double targetWindowSize = 5.5;
//
//    MsFrameScoretron msFrameScoretron;
//    e = msFrameScoretron.init(
//            PythiaParameterReader::genericPythiaParametersForTests(),
//            msDataFilePath,
//            fragLibPath,
//            fragLibBackgroundPath,
//            uniqueMsInfoScanKey,
//            {target - targetWindowSize, target + targetWindowSize}
//            );
//    QCOMPARE(e, eNoError);
//
//    QVector<ScoredCandidate> scoredCandidates;
//    e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates);
//    QCOMPARE(e, eNoError);

}

namespace {

    Err buildFragLib(
            const QString &fragLibPath,
            double targetMz,
            double isoWinMin,
            double isoWinMax,
            int topNMs2Ions,
            QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptide
            ) {

        ERR_INIT

        for (int charge : {2, 3}) {

            QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptide;

            const double massMin = BiophysicalCalcs::calculateMassFromThomson(targetMz - isoWinMin, charge, 0);
            const double massMax = BiophysicalCalcs::calculateMassFromThomson(targetMz + isoWinMin, charge, 0);

            e = FragLibReader::getMS2Ions(
                    fragLibPath,
                    AminoAcids(),
                    massMin,
                    massMax,
                    100.0,
                    2000.0,
                    topNMs2Ions,
                    true,
                    &peptideSequenceChargeKeyVsCandidatePeptide
                    ); ree;

            for (auto it = peptideSequenceChargeKeyVsCandidatePeptide.begin(); it != peptideSequenceChargeKeyVsCandidatePeptide.end(); it++) {

                const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
                const CandidatePeptide &candidatePeptide = it.value();

                PeptideStringWithMods peptideStringWithMods;
                Charge charge1;
                e = FragLibReader::peptideStringWithModsFromPeptideSequenceChargeKey(
                        peptideSequenceChargeKey,
                        &peptideStringWithMods,
                        &charge1
                        );

                const double mz = BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, charge1);

                if (mz < targetMz - isoWinMin || mz > targetMz + isoWinMax) {
                    continue;
                }

                peptideStringWithModsVsCandidatePeptide->insert(peptideStringWithMods, candidatePeptide);
            }
        }

        ERR_RETURN
    }


}//namespace
void MsFrameScoretronTests::scoreCandidatesRecalTest() {

    ERR_INIT

    const QString iRTReCalFilePath
        = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.iRT");

    const QString msDataFilePath
            = QStringLiteral("/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Desktop/Testing/LatestStuff/predicted_lib.speclib.fragLib");
//            = QStringLiteral("/home/anichols/Desktop/Libraries/2022.08.31UP000005640_9606.target.decoy.contam.human_plasma.fasta.csv.fragLib");

    const double target = 454.957;
    const UniqueMsInfoScanKey uniqueMsInfoScanKey = "454957";
    const double isoWinWidth = 5.5;

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderParquet.getMsScanInfos();
    for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {
        scanNumberVsScanTime.insert(it.key(), it.value().scanTime);
    }

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrame;
    e = msReaderParquet.collateTandemPrecursorTargetsDIA(&diaTargetFrame);
    QCOMPARE(e, eNoError);

    const int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
    e = msReaderParquet.getScanPoints(msLevel, &scanNumberVsScanPointsMS1);
    QCOMPARE(e, eNoError);

    const PythiaParameters pythiaParameters = PythiaParameterReader::genericPythiaParametersForTests();
    pythiaParameters.print();

    QMap<PeptideStringWithMods, CandidatePeptide> peptideStringWithModsVsCandidatePeptide;
    e = buildFragLib(
            fragLibPath,
            target,
            isoWinWidth,
            isoWinWidth,
            pythiaParameters.topNMs2Ions,
            &peptideStringWithModsVsCandidatePeptide
            );
    QCOMPARE(e, eNoError);

    MsFrameScoretron msFrameScoretron;
    e = msFrameScoretron.init(
            uniqueMsInfoScanKey,
            pythiaParameters,
            pythiaParameters.topNMs2Ions,
            diaTargetFrame.value(uniqueMsInfoScanKey),
            scanNumberVsScanPointsMS1,
            peptideStringWithModsVsCandidatePeptide,
            scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

    QVector<ScoredCandidate> scoredCandidates;
    e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(MsFrameScoretronTests)

#include "MsFrameScoretronTests.moc"
