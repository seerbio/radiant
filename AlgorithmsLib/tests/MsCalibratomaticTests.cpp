//
// Created by anichols on 11/07/2021.
//

#include "MsCalibratomatic.h"
#include "PythiaParameterReader.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>

class MsCalibratomaticTests : public QObject
{

Q_OBJECT

public:
    MsCalibratomaticTests();
    ~MsCalibratomaticTests() override = default;


private Q_SLOTS:

    void execTests();


private:

    PythiaParameters pythiaParameters() {

        PythiaParameters pythiaParameters;

        pythiaParameters.returnPSMTopN = 1;
        pythiaParameters.maxTandemPointCount = 2;
        pythiaParameters.ms2ExtractionWidthPPM = 12.0;
        pythiaParameters.precursorExtractionWindowThomsons = 1.0;
        pythiaParameters.chargeStateMin = 2;
        pythiaParameters.chargeStateMax = 3;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
        );

        return pythiaParameters;
    }

    Err buildPeptidesWithModsVsPeptideSequences(
            const QString &m_pepLibUri,
            QMap<PeptideStringWithMods, PeptideSequence> *peptidesWithModsVsPeptideSequences
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(m_pepLibUri); ree;

        QVector<PeptideSequence> peptideSequences;
        e = ParquetReader::read(
                m_pepLibUri,
                &peptideSequences
        );ree

        e = ErrorUtils::isNotEmpty(peptideSequences); ree;


        for (const PeptideSequence &ps : peptideSequences) {
            peptidesWithModsVsPeptideSequences->insert(ps.sequence, ps);
        }

        ERR_RETURN
    }

};

MsCalibratomaticTests::MsCalibratomaticTests() : QObject(){}

void MsCalibratomaticTests::execTests() {

    ERR_INIT

    const QString scoreVectorsFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores");

    const QString scansFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans");

    const QString pepLibFilePath
            = QDir(qApp->applicationDirPath()).filePath("2022_02_22_Homo_sapiens_UP000005640.fasta.pepLib");

    const QMap<QString, QString> filePaths = {
            {scoreVectorsFilePath, scansFilePath}
    };

    QMap<PeptideStringWithMods, PeptideSequence> peptidesWithModsVsPeptideSequences;
    e = buildPeptidesWithModsVsPeptideSequences(
            pepLibFilePath,
            &peptidesWithModsVsPeptideSequences
            );
    QCOMPARE(e, eNoError);

//    MsCalibratomatic calibratomatic;
//    e = calibratomatic.init(
//            pythiaParameters(),
//            peptidesWithModsVsPeptideSequences
//            );
//    QCOMPARE(e, eNoError);
//
//    e = calibratomatic.exec(filePaths);
//    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsCalibratomaticTests)
#include "MsCalibratomaticTests.moc"
