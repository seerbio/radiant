#include "TargetDecoyCandidatePair.h"

#include <QSemaphore>
#include <QtTest>
#include <thread>
#include <vector>

namespace {
FragLibReaderRow libraryRow() {
    FragLibReaderRow row;
    row.peptideSequenceChargeKey = "VTWLSPTNK|2";
    row.precursorCharge = 2;
    row.mass = 1045.5;
    QStringList labels;
    for (int i = 0; i < 24; ++i) {
        row.mzVals.push_back(300.0f + 20.0f * i);
        row.intensityVals.push_back(1.0f / (i + 1));
        labels.push_back(QString(i % 2 ? "b%1" : "y%1").arg(2 + i % 7));
    }
    row.ionLabels = labels.join(";");
    return row;
}

void compareIons(const QVector<MS2Ion> &actual, const QVector<MS2Ion> &expected) {
    QCOMPARE(actual.size(), expected.size());
    for (int i = 0; i < actual.size(); ++i) {
        QCOMPARE(actual.at(i).mz, expected.at(i).mz);
        QCOMPARE(actual.at(i).intensity, expected.at(i).intensity);
        QCOMPARE(actual.at(i).ionLabel, expected.at(i).ionLabel);
        QCOMPARE(actual.at(i).charge, expected.at(i).charge);
        QCOMPARE(actual.at(i).rank, expected.at(i).rank);
    }
}
}

class TargetDecoyCandidatePairConcurrencyTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void concurrentColdCacheAccess_data();
    void concurrentColdCacheAccess();
    void copiesKeepIndependentCacheInvalidation();
};

void TargetDecoyCandidatePairConcurrencyTests::concurrentColdCacheAccess_data() {
    QTest::addColumn<bool>("terminalShift");
    QTest::addColumn<bool>("sharedDecoySequence");
    QTest::newRow("penultimate") << false << false;
    QTest::newRow("terminal") << true << false;
    QTest::newRow("penultimate-collision") << false << true;
    QTest::newRow("terminal-collision") << true << true;
}

void TargetDecoyCandidatePairConcurrencyTests::concurrentColdCacheAccess() {
    QFETCH(bool, terminalShift);
    QFETCH(bool, sharedDecoySequence);
    auto row = libraryRow();
    const auto configure = [&](TargetDecoyCandidatePair *candidate) {
        candidate->setFragLibReaderRowPntr(&row);
        candidate->setDecoyFragmentShiftMode(terminalShift
            ? DecoyFragmentShiftMode::ShiftTerminalByPenultimate
            : DecoyFragmentShiftMode::ShiftPenultimate);
        candidate->decoySharesSequenceWithOtherTarget(sharedDecoySequence);
    };
    TargetDecoyCandidatePair reference(PeptideStringWithMods("VTWLSPTNK"), 0.0f);
    configure(&reference);
    const auto expectedTargets = reference.ms2IonsTarget();
    const auto expectedDecoys = reference.ms2IonsDecoy();
    constexpr int workers = 12;
    // Overlapping DIA windows can score the same candidate concurrently.
    // Start with empty caches each round and interleave target/decoy first use.
    for (int round = 0; round < 64; ++round) {
        TargetDecoyCandidatePair candidate(PeptideStringWithMods("VTWLSPTNK"), 0.0f);
        configure(&candidate);
        QSemaphore ready, start;
        std::vector<QVector<MS2Ion>> targets(workers), decoys(workers);
        std::vector<std::thread> threads;
        for (int worker = 0; worker < workers; ++worker) {
            threads.emplace_back([&, worker] {
                ready.release();
                start.acquire();
                for (int repeat = 0; repeat < 4; ++repeat) {
                    if (worker % 2) {
                        decoys[worker] = candidate.ms2IonsDecoy();
                        targets[worker] = candidate.ms2IonsTarget();
                    } else {
                        targets[worker] = candidate.ms2IonsTarget();
                        decoys[worker] = candidate.ms2IonsDecoy();
                    }
                }
            });
        }
        ready.acquire(workers);
        start.release(workers);
        for (auto &thread : threads) thread.join();
        for (int worker = 0; worker < workers; ++worker) {
            compareIons(targets[worker], expectedTargets);
            compareIons(decoys[worker], expectedDecoys);
        }
    }
}

void TargetDecoyCandidatePairConcurrencyTests::copiesKeepIndependentCacheInvalidation() {
    auto row = libraryRow();
    TargetDecoyCandidatePair original(PeptideStringWithMods("VTWLSPTNK"), 0.0f);
    original.setFragLibReaderRowPntr(&row);
    const auto expectedTargets = original.ms2IonsTarget();
    const auto expectedDecoys = original.ms2IonsDecoy();
    auto copy = original;
    copy.decoySharesSequenceWithOtherTarget(true);
    const auto shifted = copy.ms2IonsDecoy();
    QCOMPARE(shifted.size(), expectedDecoys.size());
    for (int i = 0; i < shifted.size(); ++i)
        QCOMPARE(shifted.at(i).mz, float(expectedDecoys.at(i).mz + 0.1));
    compareIons(original.ms2IonsTarget(), expectedTargets);
    compareIons(original.ms2IonsDecoy(), expectedDecoys);
    auto changedRow = row;
    changedRow.mzVals[0] += 5.0f;
    copy.setFragLibReaderRowPntr(&changedRow);
    QVERIFY(copy.ms2IonsTarget().first().mz != expectedTargets.first().mz);
    compareIons(original.ms2IonsTarget(), expectedTargets);
}

QTEST_MAIN(TargetDecoyCandidatePairConcurrencyTests)
#include "TargetDecoyCandidatePairConcurrencyTests.moc"
