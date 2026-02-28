//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "FragLibTsvReader.h"
#include "AminoAcids.h"
#include "TargetDecoyCandidatePair.h"


#include <QtTest/QtTest>

class FragLibTsvReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibTsvReaderTests() = default;
    ~FragLibTsvReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaderRowsTest();

    static void compareTest();

    static void inferIonLabelsStandardModeRoundTripTest();
    static void inferIonLabelsAlternativeModeRoundTripTest();
};

void FragLibTsvReaderTests::getFragLibReaderRowsTest() {

    ERR_INIT

    const QString &testFile = QDir(qApp->applicationDirPath()).filePath("tsv_library_trunc.tsv");

    FragLibTsvReader fragLibTsvReader;

    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(fragLibReaderRows.size(), 57);

    const FragLibReaderRow &fragLibReaderRow = fragLibReaderRows.at(47);
    qDebug() << fragLibReaderRow.peptideSequenceChargeKey
             << fragLibReaderRow.precursorCharge
             << fragLibReaderRow.mzVals
             << fragLibReaderRow.isDecoy
             << fragLibReaderRow.intensityVals
             << fragLibReaderRow.mass
             << fragLibReaderRow.iRT
             << fragLibReaderRow.ionLabels;

    QCOMPARE(fragLibReaderRow.precursorCharge, 1);
    QCOMPARE(fragLibReaderRow.mzVals.size(), 7);
    QCOMPARE(fragLibReaderRow.intensityVals.size(), 7);
    QCOMPARE(fragLibReaderRow.isDecoy, 0);
    QCOMPARE(static_cast<int>(fragLibReaderRow.mass), 600);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iRT), -38);
    QCOMPARE(fragLibReaderRow.ionLabels.size(), 20);

    QCOMPARE(e, eNoError);

}

void FragLibTsvReaderTests::compareTest() {

    QSKIP("For testing");

    ERR_INIT

    QString testFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv");
    // testFile = "/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib_trunc.tsv";
    // testFile = "/home/andrewnichols/Desktop/Data/Libraries/Test.tsv";
    // testFile = "/home/andrewnichols/Downloads/pythia-firstpass-lib.tsv";
    testFile = "/home/andrewnichols/Downloads/lib-mbr-2025-02-06.tsv";

    FragLibTsvReader fragLibTsvReader;
    QList<FragLibReaderRow> fragLibReaderTSVRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderTSVRows
    );


    // for (int i = 0; i < fragLibReaderTSVRows.size(); i++) {
    //
    //     const FragLibReaderRow &tsv = fragLibReaderTSVRows.at(i);
    //     qDebug()
    //     << i
    //     << tsv.peptideSequenceChargeKey
    //     << tsv.mzVals
    //     << tsv.intensityVals
    //     << tsv.ionLabels;
    //
    // }
}

namespace {
    constexpr int kCharge = 1;
    // first/second/middle/penultimate/last fragment indices (1,2,4,6,7) across b/y.
    const QStringList kRepresentativeIonLabels = {"y1", "b2", "y2", "b4", "y4", "b6", "y6", "b7", "y7"};
    const QString kExpectedLabels = kRepresentativeIonLabels.join(S_GLOBAL_SETTINGS.SEPARATOR);
    const PeptideStringWithMods kTargetPeptide = PeptideStringWithMods("WWTVSLPW");

    QVector<float> mzValsForIonLabels(
            const PeptideStringWithMods &peptideStringWithMods,
            const QStringList &ionLabels
            ) {
        AminoAcids aminoAcids;
        const QVector<double> bSeries = peptideStringWithMods.bSeries(kCharge, aminoAcids);
        const QVector<double> ySeries = peptideStringWithMods.ySeries(kCharge, aminoAcids);

        QVector<float> mzVals;
        mzVals.reserve(ionLabels.size());
        for (const QString &ionLabel : ionLabels) {
            const bool isB = ionLabel.startsWith('b');
            const int ionIndex = ionLabel.mid(1).toInt();
            const double mz = isB ? bSeries.at(ionIndex - 1) : ySeries.at(ionIndex - 1);
            mzVals.push_back(static_cast<float>(mz));
        }
        return mzVals;
    }

    FragLibReaderRow buildFragLibRowForIons(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<float> &mzVals,
            const QStringList &ionLabels
            ) {
        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.peptideSequenceChargeKey = peptideStringWithMods + "|" + QString::number(kCharge);
        fragLibReaderRow.precursorCharge = kCharge;
        fragLibReaderRow.isDecoy = false;
        fragLibReaderRow.ionLabels = ionLabels.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.mzVals = mzVals;
        fragLibReaderRow.intensityVals = QVector<float>(mzVals.size(), 1.0f);
        return fragLibReaderRow;
    }

    QVector<float> mzValsFromIonsByLabel(
            const QVector<MS2Ion> &ions,
            const QStringList &ionLabels
            ) {
        QMap<QString, float> ionLabelToMz;
        for (const MS2Ion &ion : ions) {
            ionLabelToMz.insert(ion.ionLabel, ion.mz);
        }

        QVector<float> mzVals;
        mzVals.reserve(ionLabels.size());
        for (const QString &ionLabel : ionLabels) {
            if (!ionLabelToMz.contains(ionLabel)) {
                return {};
            }
            mzVals.push_back(ionLabelToMz.value(ionLabel));
        }
        return mzVals;
    }

    QVector<float> decoyMzValsForMode(
            const PeptideStringWithMods &targetPeptideStringWithMods,
            const QStringList &ionLabels,
            const DecoyFragmentShiftMode shiftMode
            ) {
        const QVector<float> targetMzVals = mzValsForIonLabels(targetPeptideStringWithMods, ionLabels);
        FragLibReaderRow fragLibReaderRow = buildFragLibRowForIons(
                targetPeptideStringWithMods,
                targetMzVals,
                ionLabels
                );

        TargetDecoyCandidatePair tdcp(targetPeptideStringWithMods, 0.0f);
        tdcp.setFragLibReaderRowPntr(&fragLibReaderRow);
        tdcp.setDecoyFragmentShiftMode(shiftMode);
        return mzValsFromIonsByLabel(tdcp.ms2IonsDecoy(), ionLabels);
    }

    void compareMzVectorsNear(
            const QVector<float> &lhs,
            const QVector<float> &rhs,
            const float tolerance = 0.001f
            ) {
        QCOMPARE(lhs.size(), rhs.size());
        for (int i = 0; i < lhs.size(); ++i) {
            QVERIFY2(
                    std::abs(lhs.at(i) - rhs.at(i)) <= tolerance,
                    qPrintable(QString("m/z mismatch at index %1: %2 vs %3").arg(i).arg(lhs.at(i)).arg(rhs.at(i)))
                    );
        }
    }
}

void FragLibTsvReaderTests::inferIonLabelsStandardModeRoundTripTest() {
    const PeptideStringWithMods targetPeptideStringWithMods = kTargetPeptide;
    const PeptideStringWithMods mutatedPeptideStringWithMods
            = AminoAcids::mutatePenultimatePeptideResidues(targetPeptideStringWithMods);

    // In standard mode, emitted decoy m/z from output logic should match mutated-sequence b/y masses.
    const QVector<float> mzValsFromMutatedSeries
            = mzValsForIonLabels(mutatedPeptideStringWithMods, kRepresentativeIonLabels);
    const QVector<float> mzValsFromDecoyOutput = decoyMzValsForMode(
            targetPeptideStringWithMods,
            kRepresentativeIonLabels,
            DecoyFragmentShiftMode::ShiftPenultimate
            );
    QVERIFY2(!mzValsFromDecoyOutput.isEmpty(), "Did not find all representative ion labels in decoy ions.");
    compareMzVectorsNear(mzValsFromMutatedSeries, mzValsFromDecoyOutput);

    QString ionLabels;
    const Err e = FragLibTsvReader::inferIonLabelsForTest(
            mzValsFromDecoyOutput,
            mutatedPeptideStringWithMods,
            true,
            kCharge,
            false,
            &ionLabels
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(ionLabels, kExpectedLabels);
}

void FragLibTsvReaderTests::inferIonLabelsAlternativeModeRoundTripTest() {
    const PeptideStringWithMods targetPeptideStringWithMods = kTargetPeptide; // unmutated source written to decoy generator
    const PeptideStringWithMods mutatedPeptideStringWithMods
            = AminoAcids::mutatePenultimatePeptideResidues(targetPeptideStringWithMods);

    const QVector<float> mzVals = decoyMzValsForMode(
            targetPeptideStringWithMods,
            kRepresentativeIonLabels,
            DecoyFragmentShiftMode::ShiftTerminalByPenultimate
            );
    QVERIFY2(!mzVals.isEmpty(), "Did not find all representative ion labels in decoy ions.");

    QString ionLabels;
    Err e = FragLibTsvReader::inferIonLabelsForTest(
            mzVals,
            mutatedPeptideStringWithMods,
            true,
            kCharge,
            true,
            &ionLabels
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(ionLabels, kExpectedLabels);
}


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
