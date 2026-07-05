//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "FragLibTsvReader.h"
#include "AminoAcids.h"
#include "TargetDecoyCandidatePair.h"


#include <QFile>
#include <QTextStream>
#include <QtTest/QtTest>

class FragLibTsvReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibTsvReaderTests() = default;
    ~FragLibTsvReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaderRowsTest();
    static void getFragLibReaderRowsImsTimeAliasTest();
    static void getFragLibReaderRowsAverageExperimentalRtAliasTest();
    static void getFragLibReaderRowsApdAliasTest();
    static void getFragLibReaderRowsNoDecoyColumnTest();
    static void getFragLibReaderRowsFiltersInvalidSequencesTest();

    static void compareTest();

    static void inferIonLabelsStandardModeRoundTripTest();
    static void inferIonLabelsAlternativeModeRoundTripTest();
    static void inferIonLabelsAlternativeModeMultiplePreimageResidualSelectionTest();
    static void inferIonLabelsTargetWithNTermModificationTest();
    static void inferIonLabelsTargetWithCTermModificationTest();
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
    QVERIFY(qAbs(fragLibReaderRow.iM - 1.013057) < 0.0001);
    QCOMPARE(fragLibReaderRow.ionLabels.size(), 20);

    QCOMPARE(e, eNoError);

}

void FragLibTsvReaderTests::getFragLibReaderRowsImsTimeAliasTest() {

    ERR_INIT

    const QString testFile = QDir(qApp->applicationDirPath()).filePath("FragLibTsvReaderTests.ims_time.tsv");
    QFile file(testFile);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text));

    QTextStream out(&file);
    out << "PrecursorMz\tProductMz\tTr_recalibrated\tIMS time\tLibraryIntensity\tDecoy\tModifiedPeptide\tPrecursorCharge\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\n";
    out << "500.0\t100.0\t10.0\t1.2345\t1000.0\t0\tACDE\t2\tb\t1\t1\n";
    out << "500.0\t200.0\t10.0\t1.2345\t500.0\t0\tACDE\t2\ty\t1\t1\n";
    file.close();

    FragLibTsvReader fragLibTsvReader;
    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(e, eNoError);
    QCOMPARE(fragLibReaderRows.size(), 1);
    QVERIFY(qAbs(fragLibReaderRows.front().iM - 1.2345) < 0.0001);

    QFile::remove(testFile);
}

void FragLibTsvReaderTests::getFragLibReaderRowsAverageExperimentalRtAliasTest() {

    ERR_INIT

    const QString testFile = QDir(qApp->applicationDirPath()).filePath("FragLibTsvReaderTests.average_experimental_rt.tsv");
    QFile file(testFile);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text));

    QTextStream out(&file);
    out << "PrecursorMz\tProductMz\tProteinId\tNormalizedRetentionTime\tAverageExperimentalRetentionTime\tPrecursorIonMobility\tLibraryIntensity\tdecoy\tModifiedPeptideSequence\tPrecursorCharge\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\n";
    out << "500.0\t100.0\tP12345\t12.0\t456.7\t1.2345\t1000.0\t0\tACDEFGH\t2\tb\t1\t1\n";
    out << "500.0\t200.0\tP12345\t12.0\t456.7\t1.2345\t500.0\t0\tACDEFGH\t2\ty\t1\t1\n";
    file.close();

    FragLibTsvReader fragLibTsvReader;
    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(e, eNoError);
    QCOMPARE(fragLibReaderRows.size(), 1);
    QVERIFY(qAbs(fragLibReaderRows.front().iRT - 456.7) < 0.0001);
    QVERIFY(qAbs(fragLibReaderRows.front().iM - 1.2345) < 0.0001);
    QCOMPARE(fragLibReaderRows.front().proteinGroups, QString("P12345"));

    QFile::remove(testFile);
}

void FragLibTsvReaderTests::getFragLibReaderRowsApdAliasTest() {

    ERR_INIT

    const QString &testFile = QDir(qApp->applicationDirPath()).filePath("tsv_library_apd_alias_trunc.tsv");

    FragLibTsvReader fragLibTsvReader;

    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(fragLibReaderRows.size(), 1);

    const FragLibReaderRow &fragLibReaderRow = fragLibReaderRows.at(0);
    qDebug() << fragLibReaderRow.peptideSequenceChargeKey
             << fragLibReaderRow.precursorCharge
             << fragLibReaderRow.mzVals
             << fragLibReaderRow.isDecoy
             << fragLibReaderRow.intensityVals
             << fragLibReaderRow.mass
             << fragLibReaderRow.iRT
             << fragLibReaderRow.ionLabels;

    QCOMPARE(fragLibReaderRow.precursorCharge, 2);
    QCOMPARE(fragLibReaderRow.mzVals.size(), 2);
    QCOMPARE(fragLibReaderRow.intensityVals.size(), 2);
    QCOMPARE(fragLibReaderRow.isDecoy, 0);
    QCOMPARE(static_cast<int>(fragLibReaderRow.mass), 997);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iRT), 10);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iM), 1);
    QCOMPARE(fragLibReaderRow.ionLabels.size(), 5);

    QCOMPARE(e, eNoError);
}

void FragLibTsvReaderTests::getFragLibReaderRowsNoDecoyColumnTest() {

    ERR_INIT

    const QString &testFile = QDir(qApp->applicationDirPath()).filePath("tsv_library_no_decoy_trunc.tsv");

    FragLibTsvReader fragLibTsvReader;

    QList<FragLibReaderRow> fragLibReaderRows;
    e = fragLibTsvReader.getFragLibReaderRows(
            testFile,
            &fragLibReaderRows
            );

    QCOMPARE(fragLibReaderRows.size(), 1);

    const FragLibReaderRow &fragLibReaderRow = fragLibReaderRows.at(0);
    qDebug() << fragLibReaderRow.peptideSequenceChargeKey
             << fragLibReaderRow.precursorCharge
             << fragLibReaderRow.mzVals
             << fragLibReaderRow.isDecoy
             << fragLibReaderRow.intensityVals
             << fragLibReaderRow.mass
             << fragLibReaderRow.iRT
             << fragLibReaderRow.ionLabels;

    QCOMPARE(fragLibReaderRow.precursorCharge, 2);
    QCOMPARE(fragLibReaderRow.mzVals.size(), 2);
    QCOMPARE(fragLibReaderRow.intensityVals.size(), 2);
    QCOMPARE(fragLibReaderRow.isDecoy, 0);
    QCOMPARE(static_cast<int>(fragLibReaderRow.mass), 997);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iRT), 10);
    QCOMPARE(static_cast<int>(fragLibReaderRow.iM), 1);
    QCOMPARE(fragLibReaderRow.ionLabels.size(), 5);

    QCOMPARE(e, eNoError);
}

void FragLibTsvReaderTests::getFragLibReaderRowsFiltersInvalidSequencesTest() {
    ERR_INIT

    const QString testFile = QDir(qApp->applicationDirPath()).filePath("FragLibTsvReaderTests.invalid_sequences.tsv");
    QFile file(testFile);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text));

    QTextStream out(&file);
    out << "PrecursorMz\tProductMz\tTr_recalibrated\tIonMobility\tLibraryIntensity\tdecoy\tModifiedPeptide\tPrecursorCharge\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\n";
    out << "500.0\t100.0\t10.0\t1.0\t1000.0\t0\tACDE\t2\tb\t1\t1\n";
    out << "500.0\t200.0\t10.0\t1.0\t500.0\t0\tACDE\t2\ty\t1\t1\n";
    out << "600.0\t110.0\t20.0\t1.1\t900.0\t0\tACDX\t2\tb\t1\t1\n";
    out << "600.0\t210.0\t20.0\t1.1\t450.0\t0\tACDX\t2\ty\t1\t1\n";
    file.close();

    QList<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(testFile, false, &fragLibReaderRows);
    QCOMPARE(e, eNoError);

    QCOMPARE(fragLibReaderRows.size(), 1);
    QCOMPARE(fragLibReaderRows.front().peptideSequenceChargeKey, QString("ACDE|2"));
    QCOMPARE(static_cast<int>(fragLibReaderRows.front().iM), 1);
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

    QVector<MS2Ion> decoyIonsForMode(
            const PeptideStringWithMods &targetPeptideStringWithMods,
            const DecoyFragmentShiftMode shiftMode
            ) {
        AminoAcids aminoAcids;
        QStringList allIonLabels = targetPeptideStringWithMods.bSeriesIonLabels({});
        allIonLabels.append(targetPeptideStringWithMods.ySeriesIonLabels({}));

        const QVector<double> targetBSeries = targetPeptideStringWithMods.bSeries(kCharge, aminoAcids);
        const QVector<double> targetYSeries = targetPeptideStringWithMods.ySeries(kCharge, aminoAcids);
        QVector<float> targetMzVals;
        targetMzVals.reserve(allIonLabels.size());
        for (int i = 0; i < targetBSeries.size(); ++i) {
            targetMzVals.push_back(static_cast<float>(targetBSeries.at(i)));
        }
        for (int i = 0; i < targetYSeries.size(); ++i) {
            targetMzVals.push_back(static_cast<float>(targetYSeries.at(i)));
        }

        FragLibReaderRow fragLibReaderRow = buildFragLibRowForIons(
                targetPeptideStringWithMods,
                targetMzVals,
                allIonLabels
                );

        TargetDecoyCandidatePair tdcp(targetPeptideStringWithMods, 0.0f);
        tdcp.setFragLibReaderRowPntr(&fragLibReaderRow);
        tdcp.setDecoyFragmentShiftMode(shiftMode);
        return tdcp.ms2IonsDecoy();
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

void FragLibTsvReaderTests::inferIonLabelsAlternativeModeMultiplePreimageResidualSelectionTest() {
    const PeptideStringWithMods sourceA("WWTVSLPW");
    const PeptideStringWithMods sourceB("WGTVSLKW");
    const PeptideStringWithMods mutatedPeptideStringWithMods
            = AminoAcids::mutatePenultimatePeptideResidues(sourceA);
    QCOMPARE(mutatedPeptideStringWithMods, AminoAcids::mutatePenultimatePeptideResidues(sourceB));

    const QVector<MS2Ion> decoyIonsA = decoyIonsForMode(sourceA, DecoyFragmentShiftMode::ShiftTerminalByPenultimate);
    const QVector<MS2Ion> decoyIonsB = decoyIonsForMode(sourceB, DecoyFragmentShiftMode::ShiftTerminalByPenultimate);

    constexpr double maxDistanceForBothMatches = 0.38;
    constexpr float residualBiasTowardA = 0.01f;
    constexpr double residualEqualityTolerance = 1e-6;
    auto inferAlt = [](const QVector<float> &mzVals, const PeptideStringWithMods &peptide, QString *ionLabels) {
        return FragLibTsvReader::inferIonLabelsForTest(
                mzVals,
                peptide,
                true,
                kCharge,
                true,
                ionLabels
                );
    };
    auto ionLabelToMzMap = [](const QVector<MS2Ion> &ions) {
        QMap<QString, float> map;
        for (const MS2Ion &ion : ions) {
            map.insert(ion.ionLabel, ion.mz);
        }
        return map;
    };
    auto totalResidualFromLabels = [](const QVector<float> &mzVals, const QString &ionLabels, const QMap<QString, float> &map) {
        const QStringList labels = ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, Qt::SkipEmptyParts);
        if (labels.size() != mzVals.size()) {
            return std::numeric_limits<double>::max();
        }

        double residual = 0.0;
        for (int i = 0; i < labels.size(); ++i) {
            if (!map.contains(labels.at(i))) {
                return std::numeric_limits<double>::max();
            }
            residual += std::abs(static_cast<double>(mzVals.at(i)) - static_cast<double>(map.value(labels.at(i))));
        }
        return residual;
    };
    const QMap<QString, float> ionLabelToMzA = ionLabelToMzMap(decoyIonsA);
    const QMap<QString, float> ionLabelToMzB = ionLabelToMzMap(decoyIonsB);

    QVector<float> ambiguousObservedMzCandidates;
    ambiguousObservedMzCandidates.reserve(decoyIonsA.size() * decoyIonsB.size());
    for (const MS2Ion &ionA : decoyIonsA) {
        for (const MS2Ion &ionB : decoyIonsB) {
            const double delta = std::abs(static_cast<double>(ionA.mz) - static_cast<double>(ionB.mz));
            if (delta > maxDistanceForBothMatches) {
                continue;
            }

            ambiguousObservedMzCandidates.push_back(static_cast<float>(
                    (static_cast<double>(ionA.mz) + static_cast<double>(ionB.mz)) / 2.0
                    + (ionA.mz >= ionB.mz ? residualBiasTowardA : -residualBiasTowardA)
                    ));
        }
    }
    QVERIFY2(!ambiguousObservedMzCandidates.isEmpty(), "No ambiguous observed m/z candidates were generated.");

    bool foundCase = false;
    for (int i = 0; i < ambiguousObservedMzCandidates.size(); ++i) {
        for (int j = i; j < ambiguousObservedMzCandidates.size(); ++j) {
            const QVector<float> mzVals = {
                    ambiguousObservedMzCandidates.at(i),
                    ambiguousObservedMzCandidates.at(j)
            };
            QString ionLabelsA;
            const Err eA = inferAlt(mzVals, sourceA, &ionLabelsA);
            if (eA != eNoError) {
                continue;
            }

            QString ionLabelsB;
            const Err eB = inferAlt(mzVals, sourceB, &ionLabelsB);
            if (eB != eNoError) {
                continue;
            }

            const double residualA = totalResidualFromLabels(mzVals, ionLabelsA, ionLabelToMzA);
            const double residualB = totalResidualFromLabels(mzVals, ionLabelsB, ionLabelToMzB);
            if (std::abs(residualA - residualB) <= residualEqualityTolerance) {
                continue;
            }

            QString ionLabelsMutated;
            const Err eMut = inferAlt(mzVals, mutatedPeptideStringWithMods, &ionLabelsMutated);
            if (eMut != eNoError) {
                continue;
            }

            QCOMPARE(ionLabelsMutated, residualA < residualB ? ionLabelsA : ionLabelsB);
            foundCase = true;
            break;
        }

        if (foundCase) {
            break;
        }
    }

    QVERIFY2(foundCase, "Did not find an ambiguous multi-preimage case with deterministic residual winner.");
}

void FragLibTsvReaderTests::inferIonLabelsTargetWithNTermModificationTest() {
    constexpr float nTermDeltaMass = 42.010565f;

    const PeptideStringWithMods basePeptideStringWithMods("ACDEF");
    const PeptideStringWithMods nTermModifiedPeptideStringWithMods("(UniMod:1)ACDEF");
    const QStringList ionLabelsToMatch = {"b1", "b2", "b4", "y1", "y2", "y4"};
    const QString expectedIonLabels = ionLabelsToMatch.join(S_GLOBAL_SETTINGS.SEPARATOR);

    // Hand-shift the fragments that physically contain the peptide N-terminus so this
    // regression does not depend on the production fragment builder handling terminal mods.
    QVector<float> mzVals = mzValsForIonLabels(basePeptideStringWithMods, ionLabelsToMatch);
    for (int i = 0; i < ionLabelsToMatch.size(); ++i) {
        if (ionLabelsToMatch.at(i).startsWith('b')) {
            mzVals[i] += nTermDeltaMass;
        }
    }

    QString inferredIonLabels;
    const Err e = FragLibTsvReader::inferIonLabelsForTest(
            mzVals,
            nTermModifiedPeptideStringWithMods,
            false,
            kCharge,
            false,
            &inferredIonLabels
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(inferredIonLabels, expectedIonLabels);
}

void FragLibTsvReaderTests::inferIonLabelsTargetWithCTermModificationTest() {
    constexpr float cTermDeltaMass = 42.010565f;

    const PeptideStringWithMods basePeptideStringWithMods("ACDEF");
    const PeptideStringWithMods cTermModifiedPeptideStringWithMods("ACDEF(42.010565)");
    const QStringList ionLabelsToMatch = {"b1", "b2", "b4", "y1", "y2", "y4"};
    const QString expectedIonLabels = ionLabelsToMatch.join(S_GLOBAL_SETTINGS.SEPARATOR);

    QVector<float> mzVals = mzValsForIonLabels(basePeptideStringWithMods, ionLabelsToMatch);
    for (int i = 0; i < ionLabelsToMatch.size(); ++i) {
        if (ionLabelsToMatch.at(i).startsWith('y')) {
            mzVals[i] += cTermDeltaMass;
        }
    }

    QString inferredIonLabels;
    const Err e = FragLibTsvReader::inferIonLabelsForTest(
            mzVals,
            cTermModifiedPeptideStringWithMods,
            false,
            kCharge,
            false,
            &inferredIonLabels
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(inferredIonLabels, expectedIonLabels);
}


QTEST_MAIN(FragLibTsvReaderTests)
#include "FragLibTsvReaderTests.moc"
