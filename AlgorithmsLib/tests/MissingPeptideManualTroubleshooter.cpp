//
// Created by anichols on 11/07/2021.
//

#include "BiophysicalCalcs.h"
#include "MsReaderParquet.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "ParallelUtils.h"
#include "PythiaParameterReader.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>


struct TroubleshootOutput : public ParquetReaderInputBase {

    QVector<double> mzScan;
    QVector<double> mzScanAll;
    QVector<double> mzPred;
    QVector<double> intzScan;
    QVector<double> intzScanAll;
    QVector<double> intzPred;

    QMap<QString, QVariant> map() override {

        using namespace PSMsReaderRowNamespace;

        return {
                {"mzScan", QVariant(qVectorToQByteArray(mzScan))},
                {"mzScanAll", QVariant(qVectorToQByteArray(mzScanAll))},
                {"mzPred", QVariant(qVectorToQByteArray(mzPred))},
                {"intzScan", QVariant(qVectorToQByteArray(intzScan))},
                {"intzScanAll", QVariant(qVectorToQByteArray(intzScanAll))},
                {"intzPred", QVariant(qVectorToQByteArray(intzPred))}
        };
    }

};


class MissingPeptideManualTroubleshooter : public QObject
{

Q_OBJECT

public:
    MissingPeptideManualTroubleshooter() = default;
    ~MissingPeptideManualTroubleshooter() override = default;


private Q_SLOTS:

    void troubleshootMissingPeptide();


};

class FrameIndexScoreResultOfTarget;

void MissingPeptideManualTroubleshooter::troubleshootMissingPeptide() {

    ERR_INIT

    const QString missingPeptide = "CQXEXNFNTXQTK";
    const double scanTime = 20.5298;

    const QString fragLibFilePath = "/home/anichols/Desktop/2022_02_22_Homo_sapiens_UP000005640.fragLib";
    const QString msDataFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq";

    //TODO needs to be completely rewritten after refactor MsFrameScoretron

//    MsReaderParquet msReaderParquet;
//    e = msReaderParquet.openFile(msDataFilePath);
//    QCOMPARE(e, eNoError);
//
//    const ScanNumber scanNumber = msReaderParquet.getNearestScanNumberFromScanTime(scanTime);
//
//    MsScanInfo msScanInfo;
//    e = msReaderParquet.getMsScanInfo(scanNumber, &msScanInfo);
//    QCOMPARE(e, eNoError);
//
//    const PythiaParameters params = PythiaParameterReader::genericPythiaParametersForTests();
//    e = ErrorUtils::isTrue(params.isValid());
//    QCOMPARE(e, eNoError);
//    params.print();
//
//    const QPair<double, double> mzTargetStartStop = {
//            msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower,
//            msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper
//            };
//
//    QMap<PeptideStringWithMods, bool> fragPredsIsDecoy;
//    QMap<PeptideStringWithMods, QVector<MS2Ion>> fragPreds;
//    e = MsFrameScoretron::buildFragIonLibForTargetMz(
//            params,
//            fragLibFilePath,
//            mzTargetStartStop,
//            &fragPreds,
//            &fragPredsIsDecoy
//    );
//    QCOMPARE(e, eNoError);
//
//    e = ErrorUtils::isNotEmpty(fragPreds);
//    QCOMPARE(e, eNoError);
//
//    const bool applySmooth2D = true;
//
//    const QString uniqueMsInfoScanKey = msScanInfo.targetScanKey();
//
//    MsFrame msFrame;
//    e = MsFrameScoretron::buildMsFrame(
//            msDataFilePath,
//            uniqueMsInfoScanKey,
//            params,
//            mzTargetStartStop,
//            applySmooth2D,
//            &msFrame
//    );
//    QCOMPARE(e, eNoError);
//
//    const double peptideMass = BiophysicalCalcs::calculatePeptideMass(missingPeptide, params.aminoAcids, {});
//    const int charge = BiophysicalCalcs::calculateChargeFromSequence(missingPeptide, params.aminoAcids, msScanInfo.precursorTargetMz);
//    const int frameIndex = msFrame.frameIndexFromScanNumber(scanNumber);
//
//    qDebug() << "Theoretical Peptide Mass" << peptideMass;
//    qDebug() << "Charge" << charge;
//    qDebug() << "Theo Peptide Mz" << BiophysicalCalcs::calculateThomsonFromMass(peptideMass, charge);
//    qDebug() << "Peptide is in frag Library" << fragPreds.contains(missingPeptide);
//    qDebug() << "ScanNumber" << scanNumber;
//    qDebug() << "FrameIndex" << frameIndex;
//    qDebug() << "Target Scan Key" << uniqueMsInfoScanKey;
//
//    const ScanNumber altScanNumber = msFrame.scanNumberFromFrameIndex(252);
//    qDebug() << "Alt ScanNumber If Diff From Above" << altScanNumber;
//
//    ScanPoints scanPoints = msFrame.getScanPointsByScanNumber(scanNumber);
//
//    bool useAltScanNumber = true;
//    if (useAltScanNumber) {
//        qDebug() << "USING ALTERNATE SCAN NUMBER";
//        qDebug() << "USING ALTERNATE SCAN NUMBER";
//        qDebug() << "USING ALTERNATE SCAN NUMBER";
//        qDebug() << "USING ALTERNATE SCAN NUMBER";
//        qDebug() << "USING ALTERNATE SCAN NUMBER";
//        scanPoints = msFrame.getScanPointsByScanNumber(altScanNumber);
//    }
//
//    const QVector<MS2Ion> predMs2Ions = fragPreds.value(missingPeptide);
//
//    const QPair<QVector<double>, QVector<double>> mzIntensityUnzip = ParallelUtils::unZip(scanPoints);
//
//    const ScanPoints scanPointsAll = msFrame.getScanPointsByScanNumber(scanNumber);
//    QPair<QVector<double>, QVector<double>> unzippedScanPointsAll = ParallelUtils::unZip(scanPoints);
//
//    TroubleshootOutput output;
//    output.mzScan = mzIntensityUnzip.first;
//    output.intzScan = mzIntensityUnzip.second;
//    output.mzScanAll = unzippedScanPointsAll.first;
//    output.intzScanAll = unzippedScanPointsAll.second;
//
//    for (const MS2Ion &ms2Ion : predMs2Ions) {
//        output.mzPred.push_back(ms2Ion.x());
//        output.intzPred.push_back(ms2Ion.y());
//    }
//
//    const QString &outputFilePath
//            = QDir(qApp->applicationDirPath()).filePath("troubleshoot.prq");
//
//    const QVector<TroubleshootOutput> outputs = {output};
//    e = ParquetReader::write(outputs, outputFilePath);
//    QCOMPARE(e, eNoError);
//
//    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> peptideScoreVecs;
//    e = MsFrameScoretron::processFrameLogic(
//            {msFrame, fragPreds},
//            params,
//            &peptideScoreVecs
//            );
//    const bool inScoreVecs = peptideScoreVecs.contains(missingPeptide);
//    qDebug() << "Peptide is in Score Vecs" << inScoreVecs;
//    if (inScoreVecs) {
//        const QVector<int> foundIons = peptideScoreVecs.value(missingPeptide).foundIonsPerFrameIndexOfTargetVec;
//        qDebug() << foundIons;
//
//        const QVector<double> scores = peptideScoreVecs.value(missingPeptide).scorePerFrameIndexOfTargetVec;
//        const double indexOfTopScore = MathUtils::findMaxIndexInVector(scores);
//        qDebug() << "Best frame index score of unfound" << indexOfTopScore;
//
//        const ScanNumber scanNumberUnfound = msFrame.scanNumberFromFrameIndex(indexOfTopScore);
//        qDebug() << "Best scan number score of unfound" << scanNumberUnfound;
//
//        const ScanPoints scanPointsUnfound = msFrame.getScanPointsByScanNumber(scanNumber);
//
//    }
//
//    MsFrameScoretron::filterByFoundMzCount(
//            params.minFoundMzPeaks,
//            &peptideScoreVecs
//    );
//    qDebug() << "Peptide is in Filtered Score Vecs" << peptideScoreVecs.contains(missingPeptide);
//
//    const QString outputFilePathFrameScores = msDataFilePath + "." + msFrame.uniqueMsInfoScanKey() + ".frameScores";
//    e = MsFrameScoretron::writeFrameTargetScoreVectors(
//            peptideScoreVecs,
//            outputFilePathFrameScores
//    );
//    QCOMPARE(e, eNoError);
//
//    MsFrameScoretronProcessormatic msFrameScoretronProcessormatic;
//    msFrameScoretronProcessormatic.init(
//            fragPreds,
//            msFrame,
//            params,
//            outputFilePathFrameScores
//    );
//    QCOMPARE(e, eNoError);
//
//    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> topCansInFrameIndexVsScore;
//    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>>> topCansInFrameIndexVsDiscScore;
//    e = msFrameScoretronProcessormatic.processLogicForFrameScores(
//            &topCansInFrameIndexVsScore,
//            &topCansInFrameIndexVsDiscScore
//    );
//    QCOMPARE(e, eNoError);
//
//    const QVector<QPair<PeptideStringWithMods, Score>> &pr = topCansInFrameIndexVsScore.value(frameIndex);
//
////    for (const QPair<PeptideStringWithMods, Score> &p : pr) {
////        qDebug() << "drewholio" << p.first << p.second;
////    }
//    qDebug() << "HOW BIG" << pr.size();
//
//    ////////
//    const QVector<double> scoreVec = peptideScoreVecs.value(missingPeptide).scorePerFrameIndexOfTargetVec;
//    qDebug() << "Max score Frame Index" << MathUtils::findMaxIndexInVector(scoreVec);
//    qDebug() << peptideScoreVecs.value(missingPeptide).foundIonsPerFrameIndexOfTargetVec;

}

QTEST_MAIN(MissingPeptideManualTroubleshooter)
#include "MissingPeptideManualTroubleshooter.moc"
