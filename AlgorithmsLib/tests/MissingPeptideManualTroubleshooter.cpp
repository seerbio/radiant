//
// Created by anichols on 11/07/2021.
//

#include "BiophysicalCalcs.h"
#include "MsReaderParquet.h"
#include "MsFrameScoretron.h"
#include "ParallelUtils.h"
#include "PythiaParameterReader.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>


struct TroubleshootOutput : public ParquetReaderInputBase {

    QVector<double> mzScan;
    QVector<double> mzPred;
    QVector<double> intzScan;
    QVector<double> intzPred;

    QMap<QString, QVariant> map() override {

        using namespace PSMsReaderRowNamespace;

        return {
                {"mzScan", QVariant(qVectorToQByteArray(mzScan))},
                {"mzPred", QVariant(qVectorToQByteArray(mzPred))},
                {"intzScan", QVariant(qVectorToQByteArray(intzScan))},
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

    const QString missingPeptide = "AQXXGGVDEAWAXXQGXQSR";
    const double scanTime = 29.1513;

    const QString fragLibFilePath = "/home/anichols/Desktop/2022_02_22_Homo_sapiens_UP000005640.fragLib";
    const QString msDataFilePath = "/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq";

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    const ScanNumber scanNumber = msReaderParquet.getNearestScanNumberFromScanTime(scanTime);

    MsScanInfo msScanInfo;
    e = msReaderParquet.getMsScanInfo(scanNumber, &msScanInfo);
    QCOMPARE(e, eNoError);

    const PythiaParameters params = PythiaParameterReader::genericPythiaParametersForTests();
    e = ErrorUtils::isTrue(params.isValid());
    QCOMPARE(e, eNoError);
    params.print();

    const QPair<double, double> mzTargetStartStop = {499.479, 510.479};

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

    e = ErrorUtils::isNotEmpty(fragPreds);
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

    for (int i = 0; i < 400; i++) {
        qDebug() << i << msFrame.scanNumberFromFrameIndex(i);
    }

    const double peptideMass = BiophysicalCalcs::calculatePeptideMass(missingPeptide, params.aminoAcids, {});
    const int charge = BiophysicalCalcs::calculateChargeFromSequence(missingPeptide, params.aminoAcids, msScanInfo.precursorTargetMz);
    qDebug() << "Theoretical Peptide Mass" << peptideMass;
    qDebug() << "Charge" << charge;
    qDebug() << "Theo Peptide Mz" << BiophysicalCalcs::calculateThomsonFromMass(peptideMass, charge);
    qDebug() << "Peptide is in frag Library" << fragPreds.contains(missingPeptide); //Why is this peptide not in frag lib
    qDebug() << "ScanNumber" << scanNumber;
    qDebug() << "FrameIndex" << msFrame.frameIndexFromScanNumber(scanNumber); //Why is frame index == 0, i.e, this scan number is not in Frame even though it should be
    qDebug() << "Target Scan Key" << uniqueMsInfoScanKey;

    const ScanPoints scanPoints = msFrame.getScanPointsByScanNumber(scanNumber);
    const QVector<MS2Ion> predMs2Ions = fragPreds.value(missingPeptide);

    const QPair<QVector<double>, QVector<double>> mzIntensityUnzip = ParallelUtils::unZip(scanPoints);

    TroubleshootOutput output;
    output.mzScan = mzIntensityUnzip.first;
    output.intzScan = mzIntensityUnzip.second;

    for (const MS2Ion &ms2Ion : predMs2Ions) {
        output.mzPred.push_back(ms2Ion.mz);
        output.intzPred.push_back(ms2Ion.intensity);
    }

    const QString &outputFilePath
            = QDir(qApp->applicationDirPath()).filePath("troubleshoot.prq");

    const QVector<TroubleshootOutput> outputs = {output};
    e = ParquetReader::write(outputs, outputFilePath);
    QCOMPARE(e, eNoError);

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> peptideScoreVecs;
    e = MsFrameScoretron::processFrameLogic(
            {msFrame, fragPreds},
            params,
            &peptideScoreVecs
            );

    qDebug() << peptideScoreVecs.value(missingPeptide).foundIonsPerFrameIndexOfTargetVec;

}

QTEST_MAIN(MissingPeptideManualTroubleshooter)
#include "MissingPeptideManualTroubleshooter.moc"
