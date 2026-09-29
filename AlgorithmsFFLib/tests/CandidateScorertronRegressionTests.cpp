#include "CandidateScorertron.h"
#include "DiscriminantScoretron.h"
#include "XICPeakManager.h"

#include <QtTest/QtTest>
#include <algorithm>
#include <cmath>
#include <cstring>

class CandidateScorertronRegressionTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void rejectsPeaksWithoutLeadingFragmentSupport_data();
    void rejectsPeaksWithoutLeadingFragmentSupport();
    void preservesMs1SignalAcrossSamplingRates_data();
    void preservesMs1SignalAcrossSamplingRates();
    void preservesIndependentFragmentThresholds();
    void preservesGlobalScoresWithZeroPrefixes_data();
    void preservesGlobalScoresWithZeroPrefixes();
};

void CandidateScorertronRegressionTests::rejectsPeaksWithoutLeadingFragmentSupport_data() {
    QTest::addColumn<bool>("lowerRankedSignal");
    QTest::newRow("unsupported-leading-fragments") << true;
    QTest::newRow("no-signal") << false;
}

void CandidateScorertronRegressionTests::rejectsPeaksWithoutLeadingFragmentSupport() {
    QFETCH(bool, lowerRankedSignal);
    FragLibReaderRow library;
    library.mass = 1000.0;
    library.precursorCharge = 2;
    QStringList labels;
    for (int ion = 0; ion < 12; ++ion) {
        library.mzVals.push_back(300.0f + ion * 30.0f);
        library.intensityVals.push_back(1000.0f - ion);
        labels.push_back("p");
    }
    library.ionLabels = labels.join(S_GLOBAL_SETTINGS.SEPARATOR);
    TargetDecoyCandidatePair pair(PeptideStringWithMods("PEPTIDEK"), 0.0f);
    pair.setFragLibReaderRowPntr(&library);

    // The six lower-ranked ions form a real chromatographic peak. None of
    // the first six ions supports it, so every proposed peak is rejected
    // by the existing leading-fragment correlation check.
    QMap<ScanNumber, ScanPoints> scans;
    QMap<ScanNumber, ScanTime> times;
    for (int scan = 0; scan < 41; ++scan) {
        ScanPoints points = {ScanPoint(100.0f, 1.0f)};
        const float intensity = 1000.0f * std::exp(-std::pow((scan - 20) / 4.0f, 2.0f));
        if (lowerRankedSignal) {
            for (int ion = 6; ion < 12; ++ion)
                points.push_back(ScanPoint(library.mzVals.at(ion), intensity));
        }
        scans.insert(scan, points);
        times.insert(scan, scan * 0.02f);
    }
    QMap<ScanNumber, ScanPoints*> pointers;
    for (auto it = scans.begin(); it != scans.end(); ++it)
        pointers.insert(it.key(), &it.value());
    MsFrame frame;
    QCOMPARE(frame.init(pointers, times), eNoError);
    XICPeakManager peaks;
    QCOMPARE(peaks.init(frame, {&pair}, 12, 20.0f), eNoError);
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.maxAnchorColumnIndex = 12;
    parameters.subtractShadows = false;
    parameters.writeFullCandidateDebug = true;
    const auto features = DiscriminantScoretron::featuresOptimization();
    CandidateScorertron scorer;
    QCOMPARE(scorer.init(parameters, MsCalibratomatic(), "500000", 12, 3.9f,
                        1.0f, {{1000, {1.0f, 0.5f, 0.25f}}}, features,
                        false, &peaks, &frame, nullptr, nullptr, nullptr, nullptr), eNoError);
    CandidateScores result;
    QCOMPARE(scorer.calculateScores(pair.ms2IonsTarget(),
                                   DiscriminantScoretron::defaultWeights(features),
                                   &pair, &result), eNoError);
    QCOMPARE(result.featuresArray.at(CosineSimSum100), 0.0f);
    QCOMPARE(result.scanNumber, -1);
    QVERIFY(result.integrations.isEmpty());
    if (lowerRankedSignal)
        QVERIFY(scorer.scoringDiagnosticsSummary().contains("no_discriminant_candidate=1"));
}

void CandidateScorertronRegressionTests::preservesMs1SignalAcrossSamplingRates_data() {
    QTest::addColumn<float>("ms1Step");
    QTest::addColumn<float>("phase");
    QTest::addColumn<float>("peakShift");
    QTest::addColumn<bool>("alignByTime");
    QTest::addColumn<float>("minimumCosine");
    QTest::addColumn<float>("maximumCosine");
    QTest::newRow("slower-ms1") << .04f << 0.0f << 0.0f << false << 0.0f << 1.01f;
    QTest::newRow("same-rate") << .02f << 0.0f << 0.0f << false << 0.0f << 1.01f;
    QTest::newRow("faster-ms1") << .01f << 0.0f << 0.0f << false << 0.0f << 1.01f;
    QTest::newRow("aligned-slower") << .04f << .007f << 0.0f << true << .95f << 1.01f;
    QTest::newRow("aligned-same-rate") << .02f << .007f << 0.0f << true << .95f << 1.01f;
    QTest::newRow("aligned-faster") << .01f << .007f << 0.0f << true << .95f << 1.01f;
    QTest::newRow("shifted-slower") << .04f << .007f << .25f << true << 0.0f << .8f;
    QTest::newRow("shifted-same-rate") << .02f << .007f << .25f << true << 0.0f << .8f;
    QTest::newRow("shifted-faster") << .01f << .007f << .25f << true << 0.0f << .8f;
}

void CandidateScorertronRegressionTests::preservesMs1SignalAcrossSamplingRates() {
    QFETCH(float, ms1Step);
    QFETCH(float, phase);
    QFETCH(float, peakShift);
    QFETCH(bool, alignByTime);
    QFETCH(float, minimumCosine);
    QFETCH(float, maximumCosine);
    FragLibReaderRow library;
    library.mass = 1000.0;
    library.precursorCharge = 2;
    QStringList labels;
    for (int ion = 0; ion < 12; ++ion) {
        library.mzVals.push_back(300.0f + ion * 30.0f);
        library.intensityVals.push_back(1000.0f - ion * 40.0f);
        labels.push_back(QString(ion < 6 ? "b" : "y") + QString::number(ion % 6 + 1));
    }
    library.ionLabels = labels.join(S_GLOBAL_SETTINGS.SEPARATOR);
    TargetDecoyCandidatePair pair(PeptideStringWithMods("PEPTIDEK"), 0.0f);
    pair.setFragLibReaderRowPntr(&library);
    const auto intensityAt = [](float time) {
        return 1000.0f * std::exp(-std::pow((time - 1.0f) / .16f, 2.0f));
    };

    QMap<ScanNumber, ScanPoints> ms2Scans, ms1Scans;
    QMap<ScanNumber, ScanTime> ms2Times, ms1Times;
    for (int scan = 0; scan <= 100; ++scan) {
        const float time = scan * .02f;
        ScanPoints points = {ScanPoint(100.0f, 1.0f)};
        for (int ion = 0; ion < 12; ++ion)
            points.push_back(ScanPoint(library.mzVals.at(ion),
                                       intensityAt(time) * library.intensityVals.at(ion) / 1000.0f));
        ms2Scans.insert(scan, points);
        ms2Times.insert(scan, time);
    }
    for (int scan = 0; scan <= static_cast<int>(std::round(2.0f / ms1Step)); ++scan) {
        const float time = scan * ms1Step + phase;
        ms1Scans.insert(scan, {ScanPoint(pair.mz(false), intensityAt(time - peakShift))});
        ms1Times.insert(scan, time);
    }
    QMap<ScanNumber, ScanPoints*> ms2Pointers, ms1Pointers;
    for (auto it = ms2Scans.begin(); it != ms2Scans.end(); ++it)
        ms2Pointers.insert(it.key(), &it.value());
    for (auto it = ms1Scans.begin(); it != ms1Scans.end(); ++it)
        ms1Pointers.insert(it.key(), &it.value());
    MsFrame ms2Frame, ms1Frame;
    QCOMPARE(ms2Frame.init(ms2Pointers, ms2Times), eNoError);
    QCOMPARE(ms1Frame.init(ms1Pointers, ms1Times), eNoError);
    TurboXIC ms1Xic;
    QCOMPARE(ms1Xic.init(ms1Frame.frameIndexVsScanPoints()), eNoError);
    XICPeakManager peaks;
    QCOMPARE(peaks.init(ms2Frame, {&pair}, 12, 20.0f), eNoError);
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.subtractShadows = false;
    parameters.alignMs1ScanTimes = alignByTime;
    const auto features = DiscriminantScoretron::featuresOptimization();
    CandidateScorertron scorer;
    QCOMPARE(scorer.init(parameters, MsCalibratomatic(), "500000", 12, 3.9f,
                        2.0f, {{1000, {1.0f, 0.5f, 0.25f, 0.125f}}}, features,
                        false, &peaks, &ms2Frame, &ms1Xic, &ms1Frame, nullptr, nullptr), eNoError);
    CandidateScores result;
    QCOMPARE(scorer.calculateScores(pair.ms2IonsTarget(),
                                   DiscriminantScoretron::defaultWeights(features),
                                   &pair, &result), eNoError);
    QVERIFY(result.scanNumber >= 0);
    QVERIFY(result.featuresArray.at(Ms1IntensityFound100) > 0.0f);
    QVERIFY(result.featuresArray.at(CosineSim100MS1) > minimumCosine);
    QVERIFY(result.featuresArray.at(CosineSim100MS1) < maximumCosine);

    // Repeated candidates reuse the same scan without changing any feature.
    CandidateScores repeated;
    QCOMPARE(scorer.calculateScores(pair.ms2IonsTarget(),
                                   DiscriminantScoretron::defaultWeights(features),
                                   &pair, &repeated), eNoError);
    QCOMPARE(repeated.featuresArray.size(), result.featuresArray.size());
    QCOMPARE(std::memcmp(repeated.featuresArray.constData(), result.featuresArray.constData(),
                         result.featuresArray.size() * sizeof(float)), 0);
    QCOMPARE(repeated.scanNumber, result.scanNumber);

    // Reinitializing the same scorer must discard spectra from the prior
    // frame, even when scan numbers and underlying point addresses recur.
    for (auto it = ms2Scans.begin(); it != ms2Scans.end(); ++it) {
        it.value()[9].ry() *= .25f;
    }
    XICPeakManager changedPeaks;
    QCOMPARE(changedPeaks.init(ms2Frame, {&pair}, 12, 20.0f), eNoError);
    CandidateScorertron freshScorer;
    for (auto *item : {&scorer, &freshScorer}) {
        QCOMPARE(item->init(parameters, MsCalibratomatic(), "500000", 12, 3.9f,
                            2.0f, {{1000, {1.0f, 0.5f, 0.25f, 0.125f}}}, features,
                            false, &changedPeaks, &ms2Frame, &ms1Xic, &ms1Frame,
                            nullptr, nullptr), eNoError);
    }
    CandidateScores reinitialized, fresh;
    QCOMPARE(scorer.calculateScores(pair.ms2IonsTarget(),
                                   DiscriminantScoretron::defaultWeights(features),
                                   &pair, &reinitialized), eNoError);
    QCOMPARE(freshScorer.calculateScores(pair.ms2IonsTarget(),
                                        DiscriminantScoretron::defaultWeights(features),
                                        &pair, &fresh), eNoError);
    QCOMPARE(reinitialized.featuresArray.size(), fresh.featuresArray.size());
    QCOMPARE(std::memcmp(reinitialized.featuresArray.constData(), fresh.featuresArray.constData(),
                         fresh.featuresArray.size() * sizeof(float)), 0);
    QVERIFY(std::memcmp(reinitialized.featuresArray.constData(), result.featuresArray.constData(),
                        result.featuresArray.size() * sizeof(float)) != 0);
    QCOMPARE(reinitialized.scanNumber, fresh.scanNumber);
}

void CandidateScorertronRegressionTests::preservesIndependentFragmentThresholds() {
    FragLibReaderRow library;
    library.mass = 1000.0;
    library.precursorCharge = 2;
    QStringList labels;
    for (int ion = 0; ion < 12; ++ion) {
        library.mzVals.push_back(300.0f + ion * 30.0f);
        library.intensityVals.push_back(1000.0f - ion * 40.0f);
        labels.push_back(QString(ion < 6 ? "b" : "y") + QString::number(ion % 6 + 1));
    }
    library.ionLabels = labels.join(S_GLOBAL_SETTINGS.SEPARATOR);
    TargetDecoyCandidatePair pair(PeptideStringWithMods("PEPTIDEK"), 0.0f);
    pair.setFragLibReaderRowPntr(&library);
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.subtractShadows = false;
    parameters.maxAnchorColumnIndex = 12;
    const auto features = DiscriminantScoretron::featuresOptimization();
    const auto weights = DiscriminantScoretron::defaultWeights(features);

    for (int strongIons : {0, 3, 4, 12}) {
        QMap<ScanNumber, ScanPoints> scans;
        QMap<ScanNumber, ScanTime> times;
        for (int scan = 0; scan <= 100; ++scan) {
            const float peak = std::exp(-std::pow((scan - 50) / 8.0f, 2.0f));
            ScanPoints points = {ScanPoint(100.0f, 1.0f)};
            const QVector<int> order = {0, 1, 6, 7, 2, 3, 4, 5, 8, 9, 10, 11};
            for (int rank = 0; rank < 12; ++rank) {
                const int ion = order.at(rank);
                const float scale = rank < strongIons ? library.intensityVals.at(ion) : .05f;
                points.push_back(ScanPoint(library.mzVals.at(ion), peak * scale));
            }
            scans.insert(scan, points);
            times.insert(scan, scan * .02f);
        }
        QMap<ScanNumber, ScanPoints*> pointers;
        for (auto it = scans.begin(); it != scans.end(); ++it) pointers.insert(it.key(), &it.value());
        MsFrame frame;
        QCOMPARE(frame.init(pointers, times), eNoError);
        XICPeakManager peaks;
        QCOMPARE(peaks.init(frame, {&pair}, 12, 20.0f), eNoError);
        const auto initialize = [&](CandidateScorertron *scorer, float threshold) {
            return scorer->init(parameters, MsCalibratomatic(), "500000", 12, threshold,
                                2.0f, {{1000, {1.0f, .5f, .25f, .125f}}}, features,
                                false, &peaks, &frame, nullptr, nullptr, nullptr, nullptr);
        };
        CandidateScorertron joint;
        QCOMPARE(initialize(&joint, 3.9f), eNoError);
        for (const QVector<float> thresholds : {
                QVector<float>{3.9f, 2.9f, 4.9f}, QVector<float>{2.9f, 4.9f, 3.9f}}) {
            QVector<CandidateScores> together(thresholds.size());
            QCOMPARE(joint.calculateScoresForFragmentThresholds(
                pair.ms2IonsTarget(), weights, &pair, thresholds, &together), eNoError);
            for (int view = 0; view < thresholds.size(); ++view) {
                CandidateScorertron independent;
                QCOMPARE(initialize(&independent, thresholds.at(view)), eNoError);
                CandidateScores expected;
                QCOMPARE(independent.calculateScores(pair.ms2IonsTarget(), weights, &pair, &expected), eNoError);
                QCOMPARE(together.at(view).featuresArray.size(), expected.featuresArray.size());
                QCOMPARE(std::memcmp(together.at(view).featuresArray.constData(), expected.featuresArray.constData(),
                                     expected.featuresArray.size() * sizeof(float)), 0);
                QCOMPARE(together.at(view).scanNumber, expected.scanNumber);
                QCOMPARE(together.at(view).scanNumberStart, expected.scanNumberStart);
                QCOMPARE(together.at(view).scanNumberEnd, expected.scanNumberEnd);
            }
            if (strongIons == 3) {
                const int three = thresholds.indexOf(2.9f);
                const int four = thresholds.indexOf(3.9f);
                QVERIFY(together.at(three).scanNumber >= 0);
                QCOMPARE(together.at(four).scanNumber, -1);
            }
        }
    }
}

void CandidateScorertronRegressionTests::preservesGlobalScoresWithZeroPrefixes_data() {
    QTest::addColumn<int>("firstSignal");
    QTest::addColumn<int>("pattern");
    QTest::addColumn<int>("kernel");
    QTest::addColumn<bool>("subtractShadows");
    for (int first : {0, 33, 257, 1025}) {
        for (int pattern : {0, 1, 2}) {
            for (int kernel : {3, 5, 7}) {
                for (bool shadows : {false, true}) {
                    const QByteArray name = QString("prefix-%1-pattern-%2-kernel-%3-shadows-%4")
                        .arg(first).arg(pattern).arg(kernel).arg(shadows).toLatin1();
                    QTest::newRow(name.constData()) << first << pattern << kernel << shadows;
                }
            }
        }
    }
}

void CandidateScorertronRegressionTests::preservesGlobalScoresWithZeroPrefixes() {
    QFETCH(int, firstSignal);
    QFETCH(int, pattern);
    QFETCH(int, kernel);
    QFETCH(bool, subtractShadows);
    FragLibReaderRow library;
    library.mass = 1000.0;
    library.precursorCharge = 2;
    QStringList labels;
    for (int ion = 0; ion < 12; ++ion) {
        library.mzVals.push_back(300.0f + ion * 30.0f);
        library.intensityVals.push_back(1000.0f - ion * 40.0f);
        labels.push_back(QString(ion < 6 ? "b" : "y") + QString::number(ion % 6 + 1));
    }
    library.ionLabels = labels.join(S_GLOBAL_SETTINGS.SEPARATOR);
    TargetDecoyCandidatePair pair(PeptideStringWithMods("PEPTIDEK"), 0.0f);
    pair.setFragLibReaderRowPntr(&library);
    const auto ions = pair.ms2IonsTarget();
    const int firstPrimary = firstSignal + (pattern == 2 ? 20 : 0);
    const int width = pattern == 0 ? 12 : 41;
    const int lastFrame = firstPrimary + 4 * width;
    const auto peakAt = [=](int frame) {
        if (frame < firstPrimary || frame > firstPrimary + 3 * width) return 0.0f;
        const int local = frame - firstPrimary;
        const float first = std::exp(-std::pow((local - width) / (width / 4.0f), 2.0f));
        // A second equal-height peak exercises the full-vector tie ordering.
        const float second = pattern == 1
            ? std::exp(-std::pow((local - 2 * width) / (width / 4.0f), 2.0f)) : 0.0f;
        return 1000.0f * (first + second);
    };
    QMap<ScanNumber, ScanPoints> scans, ms1Scans;
    QMap<ScanNumber, ScanTime> times;
    for (int frame = 0; frame <= lastFrame; ++frame) {
        ScanPoints points = {ScanPoint(100.0f, 1.0f)};
        for (int ion = 0; ion < ions.size(); ++ion) {
            if (frame >= firstPrimary && frame <= firstPrimary + 3 * width) {
                const float scale = pattern == 2 && ion >= 3 ? .00001f
                                    : library.intensityVals.at(ion) / 1000.0f;
                points.push_back(ScanPoint(ions.at(ion).mz, peakAt(frame) * scale));
            }
            // Shadow evidence starts before the primary evidence in pattern 2.
            if (frame >= firstSignal && frame <= firstPrimary + 3 * width) {
                const float shadow = pattern == 2 ? .000001f : peakAt(frame) * .015f;
                points.push_back(ScanPoint(
                    ions.at(ion).mz - S_GLOBAL_SETTINGS.ISO_DIFF / ions.at(ion).charge, shadow));
            }
        }
        std::sort(points.begin(), points.end(),
                  [](const ScanPoint &a, const ScanPoint &b) { return a.x() < b.x(); });
        const ScanNumber scan = frame * 3 + 17;
        scans.insert(scan, points);
        ms1Scans.insert(scan, {ScanPoint(pair.mz(false), peakAt(frame))});
        times.insert(scan, frame * .02f);
    }
    QMap<ScanNumber, ScanPoints*> pointers, ms1Pointers;
    for (auto it = scans.begin(); it != scans.end(); ++it) pointers.insert(it.key(), &it.value());
    for (auto it = ms1Scans.begin(); it != ms1Scans.end(); ++it) ms1Pointers.insert(it.key(), &it.value());
    MsFrame frame, ms1Frame;
    QCOMPARE(frame.init(pointers, times), eNoError);
    QCOMPARE(ms1Frame.init(ms1Pointers, times), eNoError);
    TurboXIC ms1Xic;
    QCOMPARE(ms1Xic.init(ms1Frame.frameIndexVsScanPoints()), eNoError);
    XICPeakManager peaks;
    QCOMPARE(peaks.init(frame, {&pair}, 12, 20.0f), eNoError);
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.filterLengthMS2 = kernel;
    parameters.subtractShadows = subtractShadows;
    const auto features = DiscriminantScoretron::featuresOptimization();
    const auto weights = DiscriminantScoretron::defaultWeights(features);

    const auto compare = [](const CandidateScores &actual, const CandidateScores &expected) {
        QCOMPARE(actual.featuresArray.size(), expected.featuresArray.size());
        QCOMPARE(std::memcmp(actual.featuresArray.constData(), expected.featuresArray.constData(),
                             expected.featuresArray.size() * sizeof(float)), 0);
        QCOMPARE(actual.integrations.size(), expected.integrations.size());
        if (!expected.integrations.isEmpty()) {
            QCOMPARE(std::memcmp(actual.integrations.constData(), expected.integrations.constData(),
                                 expected.integrations.size() * sizeof(float)), 0);
        }
        QCOMPARE(actual.frameIndex, expected.frameIndex);
        QCOMPARE(actual.frameIndexStart, expected.frameIndexStart);
        QCOMPARE(actual.frameIndexEnd, expected.frameIndexEnd);
        QCOMPARE(actual.scanNumber, expected.scanNumber);
        QCOMPARE(actual.scanNumberStart, expected.scanNumberStart);
        QCOMPARE(actual.scanNumberEnd, expected.scanNumberEnd);
        for (auto field : {&CandidateScores::scanTime, &CandidateScores::scanTimeStart,
                           &CandidateScores::scanTimeEnd, &CandidateScores::scanTimePredicted,
                           &CandidateScores::empiricalIRT, &CandidateScores::imDriftTime}) {
            QCOMPARE(std::memcmp(&(actual.*field), &(expected.*field), sizeof(float)), 0);
        }
        for (auto field : {&CandidateScores::classifierScore, &CandidateScores::discriminantScore,
                           &CandidateScores::qValue, &CandidateScores::decoyRatio,
                           &CandidateScores::precursorQValue, &CandidateScores::peptideQValue,
                           &CandidateScores::proteinQValue}) {
            QCOMPARE(std::memcmp(&(actual.*field), &(expected.*field), sizeof(double)), 0);
        }
        QCOMPARE(actual.ionMobilityIndex, expected.ionMobilityIndex);
        QCOMPARE(actual.ionMobilityIndexStart, expected.ionMobilityIndexStart);
        QCOMPARE(actual.ionMobilityIndexEnd, expected.ionMobilityIndexEnd);
        QCOMPARE(actual.ionLabels, expected.ionLabels);
    };
    for (float stop : {0.0f, .2f}) {
        parameters.stopThresholdFractionMS2 = stop;
        CandidateScorertron compact, original;
        const auto initialize = [&](CandidateScorertron *scorer, bool debug) {
            auto settings = parameters;
            // The diagnostic path retains the original global matrices.
            settings.writeFullCandidateDebug = debug;
            return scorer->init(settings, MsCalibratomatic(), "500000", 12, 3.9f,
                                lastFrame * .02f, {{1000, {1.0f, .5f, .25f, .125f}}},
                                features, false, &peaks, &frame, &ms1Xic, &ms1Frame,
                                nullptr, nullptr);
        };
        QCOMPARE(initialize(&compact, false), eNoError);
        QCOMPARE(initialize(&original, true), eNoError);
        for (const QVector<float> thresholds : {QVector<float>{3.9f, 2.9f},
                                                QVector<float>{2.9f, 3.9f}}) {
            QVector<CandidateScores> actual(2), expected(2);
            QCOMPARE(compact.calculateScoresForFragmentThresholds(
                ions, weights, &pair, thresholds, &actual), eNoError);
            QCOMPARE(original.calculateScoresForFragmentThresholds(
                ions, weights, &pair, thresholds, &expected), eNoError);
            for (int view = 0; view < 2; ++view) {
                compare(actual.at(view), expected.at(view));
                if (pattern != 2 || thresholds.at(view) == 2.9f) {
                    QVERIFY(expected.at(view).frameIndex >= firstPrimary);
                    QCOMPARE(expected.at(view).scanNumber, expected.at(view).frameIndex * 3 + 17);
                }
            }
        }
    }
}

QTEST_MAIN(CandidateScorertronRegressionTests)
#include "CandidateScorertronRegressionTests.moc"
