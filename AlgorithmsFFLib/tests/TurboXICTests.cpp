//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "MsReaderMzML.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtTest/QtTest>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>

class TurboXICTests : public QObject
{
    Q_OBJECT

public:
    TurboXICTests() = default;
    ~TurboXICTests() override = default;

private Q_SLOTS:

    void initTest();
    void extractPointsTest();
    void scanRestrictedQueryPreservesOrderAndValues();
    void indexedQueriesMatchOriginalTree();
    void massDirectoryPreservesBoundariesAndFallbacks();
    void turboXICUtility();


private:

    QMap<int, QVector<PointFF>> buildPoints() const;

};


QMap<int, QVector<PointFF>> TurboXICTests::buildPoints() const {

    QMap<int, QVector<PointFF>> points;
    points.insert(1, {PointFF(100.11, 100.1), PointFF(200.11, 200.1)});
    points.insert(2, {PointFF(100.12, 100.2), PointFF(200.12, 200.2)});
    points.insert(3, {PointFF(100.13, 100.3), PointFF(200.13, 200.3)});
    points.insert(4, {PointFF(100.14, 100.4), PointFF(200.14, 200.4)});
    points.insert(5, {PointFF(100.15, 100.5), PointFF(100.151, 100.5), PointFF(200.15, 200.5)});


    return points;
}

void TurboXICTests::initTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanPoints*> emptyPoints;
    e = turboXIC.init(emptyPoints);
    QCOMPARE(e, eEmptyContainerError);

}

void TurboXICTests::extractPointsTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    //TODO fix test

    const XICPoints xicPoints = turboXIC.extractPointsXIC(100.0, 100.13);

    QCOMPARE(xicPoints.size(), 3);
    QCOMPARE(xicPoints.front().scanNumber, 1);
    QVERIFY(MathUtils::tSame(xicPoints.front().intensity, 100.1f));
    QVERIFY(MathUtils::tSame(xicPoints.back().intensity, 100.3f));

}

void TurboXICTests::turboXICUtility() {

    ERR_INIT

    QSKIP("uncomment for troubleshooting");

//    const QString msDataFilePath
//            = QStringLiteral("/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");
//
//    const MzTargetKey uniqueMsInfoScanKey = "454957";
//    double target = 454.957;
//    double window = 5.5;
//
//    MsFrame msFrame;
//    e = MsFrame::buildMsFrame(
//            msDataFilePath,
//            uniqueMsInfoScanKey,
//            {target - window, target + window},
//            &msFrame
//    );
//    QCOMPARE(e, eNoError);
//
//    e = msFrame.smoothFrame(
//            3,
//            1.0,
//            2,
//            1500.0
//            );
//    QCOMPARE(e, eNoError);
//
//    const QMap<int, QVector<PointFF>> &points = msFrame.scanNumberVsScanPoints();
//
//    TurboXIC turboXIC;
//    e = turboXIC.init(points);
//    QCOMPARE(e, eNoError);
//
//    const double mzCenter = 523.26232347;
//    const double ppmTol = 50;
//    const double massTol = MathUtils::calculatePPM(mzCenter, ppmTol);
//
//    const XICPoints xicPoints = turboXIC.extractPointsXIC(
//            mzCenter - massTol,
//            mzCenter + massTol,
//            0,
//            26000
//            );
//
//    qDebug() << xicPoints.scanNumbersVsIntensityVals.size();
//    const QVector<PointFF> vec = ParallelUtils::convertMapToPoints(xicPoints.scanNumbersVsIntensityVals);
//
//    for (const PointFF &p : vec) {
//        qDebug() << p;
//    }
//
//    e = MsUtils::writePointsToCSV(vec, "xic.csv");
//    QCOMPARE(e, eNoError);

}

void TurboXICTests::scanRestrictedQueryPreservesOrderAndValues() {
    QMap<ScanNumber, ScanPoints> storage;
    for (int scan = -2; scan < 80; ++scan) {
        ScanPoints points;
        for (int index = 0; index < 129; ++index) {
            // Repeated m/z values and multiple observations per scan make
            // any reordering of the spatial query visible.
            points.push_back(ScanPoint(100.0f + float(index % 33) * .125f,
                                      float((scan + 3) * 129 + index)));
        }
        storage.insert(scan, points);
    }
    QMap<ScanNumber, ScanPoints*> pointers;
    for (auto it = storage.begin(); it != storage.end(); ++it)
        pointers.insert(it.key(), &it.value());
    TurboXIC index;
    QCOMPARE(index.init(pointers), eNoError);
    const QVector<QPair<float, float>> massRanges{
        {99, 105}, {100, 100}, {100.125f, 102.5f}, {104, 104}, {105, 106}};
    const QVector<QPair<int, int>> scanRanges{
        {-10, 100}, {-1, 1}, {0, 0}, {3, 8}, {78, 100}, {101, 102}, {20, 10}};
    for (const auto &mass : massRanges) {
        const auto all = index.extractPointsXIC(mass.first, mass.second);
        for (const auto &scans : scanRanges) {
            for (bool inclusive : {false, true}) {
                auto expected = all;
                expected.erase(std::remove_if(expected.begin(), expected.end(),
                    [&](const XICPoint &point) {
                        return inclusive
                            ? !(scans.first <= point.scanNumber && point.scanNumber <= scans.second)
                            : !(scans.first < point.scanNumber && point.scanNumber < scans.second);
                    }), expected.end());
                const auto actual = index.extractPointsXIC(
                    mass.first, mass.second, scans.first, scans.second, inclusive);
                QCOMPARE(actual.size(), expected.size());
                for (size_t row = 0; row < actual.size(); ++row) {
                    QCOMPARE(actual[row].scanNumber, expected[row].scanNumber);
                    QCOMPARE(actual[row].ionMobilityIndex, expected[row].ionMobilityIndex);
                    QVERIFY(std::memcmp(&actual[row].mz, &expected[row].mz, sizeof(float)) == 0);
                    QVERIFY(std::memcmp(&actual[row].intensity, &expected[row].intensity,
                                        sizeof(float)) == 0);
                }
            }
        }
    }
}

void TurboXICTests::indexedQueriesMatchOriginalTree() {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;
    using Coordinate = bg::model::point<float, 1, bg::cs::cartesian>;
    using TreePoint = std::pair<Coordinate, std::pair<float, float>>;
    using Tree = bgi::rtree<TreePoint, bgi::dynamic_quadratic>;

    std::mt19937 random(666);
    QMap<ScanNumber, ScanPoints> storage;
    const QVector<ScanNumber> scans{-3, 0, 1, 7, 100, 16777217, 16777219};
    for (auto scan : scans) {
        ScanPoints points;
        for (int index = 0; index < 3000; ++index) {
            const float mz = index % 3
                ? 100.0f + float(random() % 100000) / 128.0f
                : 100.0f + float(index % 17) / 8.0f;
            float intensity = float(random());
            if (index == 0) intensity = -0.0f;
            if (index == 1) intensity = std::numeric_limits<float>::infinity();
            if (index == 2) {
                const std::uint32_t payload = 0x7fc00123u;
                std::memcpy(&intensity, &payload, sizeof(float));
            }
            points.push_back(ScanPoint(mz, intensity));
        }
        storage.insert(scan, points);
    }
    QMap<ScanNumber, ScanPoints*> pointers;
    std::vector<TreePoint> cloud;
    for (auto it = storage.begin(); it != storage.end(); ++it) {
        pointers.insert(it.key(), &it.value());
        for (const auto &point : it.value())
            cloud.emplace_back(Coordinate(point.x()),
                std::make_pair(static_cast<float>(it.key()), point.y()));
    }
    // Independent oracle: the original constructor and spatial query, before
    // either the bounded-query or the sorted-mass optimization.
    std::sort(cloud.begin(), cloud.end(), [](const TreePoint &left, const TreePoint &right) {
        return left.first.get<0>() < right.first.get<0>();
    });
    const Tree reference(cloud, bgi::dynamic_quadratic(16));
    TurboXIC indexed;
    QCOMPARE(indexed.init(pointers), eNoError);
    TurboXIC pointerInitialized;
    QCOMPARE(pointerInitialized.init(&pointers), eNoError);

    const float infinity = std::numeric_limits<float>::infinity();
    QVector<QPair<float, float>> masses{
        {-infinity, infinity}, {0, 2000}, {100, 100},
        {std::nextafter(100.0f, -infinity), std::nextafter(100.0f, infinity)},
        {std::nextafter(100.0f, infinity), 100.125f},
        {100.125f, std::nextafter(100.25f, -infinity)},
        {99, 99.5f}, {2000, 2100}};
    for (int index = 0; index < 256; ++index) {
        const float mass = cloud[random() % cloud.size()].first.get<0>();
        masses.push_back({mass, mass + float(random() % 256) / 128.0f});
    }
    const QVector<QPair<ScanNumber, ScanNumber>> ranges{
        {-3, 100}, {0, 0}, {1, 7}, {101, 102}, {100, -3},
        {16777216, 16777218}, {16777218, 16777220}};
    const auto same = [](const XICPoints &actual, const std::vector<TreePoint> &expected) {
        if (actual.size() != expected.size()) return false;
        for (size_t row = 0; row < actual.size(); ++row) {
            const float mass = expected[row].first.get<0>();
            if (actual[row].scanNumber != static_cast<ScanNumber>(expected[row].second.first)
                || actual[row].ionMobilityIndex != -1
                || std::memcmp(&actual[row].mz, &mass, sizeof(float))
                || std::memcmp(&actual[row].intensity, &expected[row].second.second, sizeof(float))) {
                return false;
            }
        }
        return true;
    };
    for (const auto &mass : masses) {
        std::vector<TreePoint> original;
        reference.query(bgi::intersects(
            bg::model::box<Coordinate>(Coordinate(mass.first), Coordinate(mass.second))),
            std::back_inserter(original));
        QVERIFY(same(indexed.extractPointsXIC(mass.first, mass.second), original));
        QVERIFY(same(pointerInitialized.extractPointsXIC(mass.first, mass.second), original));
        for (const auto &range : ranges) {
            for (const bool inclusive : {false, true}) {
                auto expected = original;
                expected.erase(std::remove_if(expected.begin(), expected.end(),
                    [&](const TreePoint &point) {
                        const auto scan = static_cast<ScanNumber>(point.second.first);
                        return inclusive
                            ? !(range.first <= scan && scan <= range.second)
                            : !(range.first < scan && scan < range.second);
                    }), expected.end());
                QVERIFY(same(indexed.extractPointsXIC(
                    mass.first, mass.second, range.first, range.second, inclusive), expected));
            }
        }
    }

    // Reinitialization must invalidate the previous mass index.
    auto smallStorage = buildPoints();
    QMap<ScanNumber, ScanPoints*> smallPointers;
    for (auto it = smallStorage.begin(); it != smallStorage.end(); ++it)
        smallPointers.insert(it.key(), &it.value());
    QCOMPARE(indexed.init(smallPointers), eNoError);
    const auto after = indexed.extractPointsXIC(100, 101);
    QCOMPARE(after.size(), size_t(6));
    QCOMPARE(after.front().scanNumber, 1);
    QCOMPARE(after.back().scanNumber, 5);
}



void TurboXICTests::massDirectoryPreservesBoundariesAndFallbacks() {
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;
    using Coordinate = bg::model::point<float, 1, bg::cs::cartesian>;
    using TreePoint = std::pair<Coordinate, std::pair<float, float>>;
    using Tree = bgi::rtree<TreePoint, bgi::dynamic_quadratic>;
    const float lowest = std::numeric_limits<float>::lowest();
    const float highest = std::numeric_limits<float>::max();
    const QVector<QVector<float>> massSets{
        {-0.25f, std::nextafter(-.125f, lowest), -.125f, std::nextafter(-.125f, highest),
         -0.0f, 0.0f, std::nextafter(.125f, lowest), .125f, std::nextafter(.125f, highest), .25f},
        {100.0f, std::nextafter(100.125f, lowest), 100.125f,
         std::nextafter(100.125f, highest), 100.5f},
        {0.0f, 1000000.0f, 1000001.0f}, // Directory size guard: retain full mass search.
        {1.0e20f, 1.0e20f},             // A single bin at a large absolute coordinate.
        {100.0f, 100.125f, 100.25f}     // Reinitialization after fallback cases.
    };
    TurboXIC index;
    for (const auto &masses : massSets) {
        QMap<ScanNumber, ScanPoints> storage;
        for (int scan = 0; scan < 7; ++scan) {
            ScanPoints points;
            for (int repeat = 0; repeat < 5; ++repeat)
                for (int i = 0; i < masses.size(); ++i)
                    points.push_back(ScanPoint(masses.at(i), float(scan * 1000 + repeat * 100 + i)));
            storage.insert(scan, points);
        }
        QMap<ScanNumber, ScanPoints*> pointers;
        std::vector<TreePoint> cloud;
        for (auto it = storage.begin(); it != storage.end(); ++it) {
            pointers.insert(it.key(), &it.value());
            for (const auto &point : it.value())
                cloud.emplace_back(Coordinate(point.x()), std::make_pair(float(it.key()), point.y()));
        }
        std::sort(cloud.begin(), cloud.end(),
                  [](const TreePoint &a, const TreePoint &b) { return a.first.get<0>() < b.first.get<0>(); });
        const Tree original(cloud, bgi::dynamic_quadratic(16));
        QCOMPARE(index.init(pointers), eNoError);
        QVector<QPair<float, float>> ranges{{lowest, highest}, {-1.0e10f, 1.0e10f}, {-0.0f, 0.0f}};
        for (float mass : masses) {
            ranges.push_back({mass, mass});
            ranges.push_back({std::nextafter(mass, lowest), mass});
            ranges.push_back({mass, std::nextafter(mass, highest)});
        }
        for (const auto &range : ranges) {
            std::vector<TreePoint> originalPoints;
            const bg::model::box<Coordinate> box(Coordinate(range.first), Coordinate(range.second));
            original.query(bgi::intersects(box), std::back_inserter(originalPoints));
            for (bool inclusive : {false, true}) {
                const auto actual = index.extractPointsXIC(range.first, range.second, 1, 5, inclusive);
                std::vector<TreePoint> expected;
                for (const auto &point : originalPoints) {
                    const int scan = static_cast<int>(point.second.first);
                    if (inclusive ? (1 <= scan && scan <= 5) : (1 < scan && scan < 5)) expected.push_back(point);
                }
                QCOMPARE(actual.size(), expected.size());
                for (std::size_t row = 0; row < expected.size(); ++row) {
                    const float mass = expected[row].first.get<0>();
                    QCOMPARE(actual[row].scanNumber, static_cast<int>(expected[row].second.first));
                    QCOMPARE(actual[row].ionMobilityIndex, -1);
                    QVERIFY(std::memcmp(&actual[row].mz, &mass, sizeof(float)) == 0);
                    QVERIFY(std::memcmp(&actual[row].intensity, &expected[row].second.second, sizeof(float)) == 0);
                }
            }
        }
    }
}

QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
