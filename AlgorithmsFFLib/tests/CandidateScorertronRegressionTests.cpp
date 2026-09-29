#include "CandidateScorertron.h"
#include "DiscriminantScoretron.h"
#include "XICPeakManager.h"

#include <QtTest/QtTest>
#include <cmath>
#include <cstring>

class CandidateScorertronRegressionTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void rejectsPeaksWithoutLeadingFragmentSupport_data();
    void rejectsPeaksWithoutLeadingFragmentSupport();
    void preservesMs1SignalAcrossSamplingRates_data();
    void preservesMs1SignalAcrossSamplingRates();
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

QTEST_MAIN(CandidateScorertronRegressionTests)
#include "CandidateScorertronRegressionTests.moc"
