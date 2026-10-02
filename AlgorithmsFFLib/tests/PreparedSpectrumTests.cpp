#include "MsUtils.h"
#include <QtTest/QtTest>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>

namespace {
// Frozen b39e6fab implementation: independent of the prepared-spectrum path.
ExtractPoints originalExtraction(
        const QVector<QPointF> &_points,
        const QVector<QPointF> &_pointsToExtract,
        double extractionPPM
        ) {

    ExtractPoints extractPointsOutput;
    extractPointsOutput.mzFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});
    extractPointsOutput.intensityFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});

    QVector<QPointF> pointsToExtract = _pointsToExtract;
    std::sort(pointsToExtract.begin(), pointsToExtract.end(), [](const QPointF &l, const QPointF &r){return l.x() < r.x();});

    QVector<QPointF> points = _points;
    std::sort(points.begin(), points.end(), [](const QPointF &l, const QPointF &r){return l.x() < r.x();});

    int currentExtractionIndex = 0;
    double extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
    double extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

    extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;

    double xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
    double xLo = extractionPointX - xPPM;
    double xHi = extractionPointX + xPPM;

    for (const QPointF &xPoint : points) {

        const double xVal = xPoint.x();

        while (xVal > xHi) {

            currentExtractionIndex++;

            if (currentExtractionIndex > _pointsToExtract.size() - 1) {
                break;
            }

            extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
            extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

            xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
            xLo = extractionPointX - xPPM;
            xHi = extractionPointX + xPPM;

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;
        }

        if (xLo <= xVal && xVal <= xHi) {

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex]
                    =  xPoint.y() > extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    ? QPointF(xVal , extractionPointX)
                    : extractPointsOutput.mzFoundVsSearched[currentExtractionIndex];

            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].rx() =  std::max(
                    xPoint.y(),
                    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    );
        }
    }

    while (extractPointsOutput.intensityFoundVsSearched.back().y() < 0) {
        extractPointsOutput.mzFoundVsSearched.pop_back();
        extractPointsOutput.intensityFoundVsSearched.pop_back();
    }

    return extractPointsOutput;
}


std::uint64_t bits(double value) {
    std::uint64_t result;
    std::memcpy(&result, &value, sizeof(result));
    return result;
}
void comparePoints(const QVector<QPointF> &actual, const QVector<QPointF> &expected) {
    QCOMPARE(actual.size(), expected.size());
    for (int i = 0; i < actual.size(); ++i) {
        QCOMPARE(bits(actual.at(i).x()), bits(expected.at(i).x()));
        QCOMPARE(bits(actual.at(i).y()), bits(expected.at(i).y()));
    }
}
}

class PreparedSpectrumTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void preservesOriginalExtraction_data() {
        QTest::addColumn<double>("ppm");
        QTest::newRow("exact-mass") << 0.0;
        QTest::newRow("narrow") << 20.0;
        QTest::newRow("overlapping-windows") << 100000.0;
    }
    void preservesOriginalExtraction() {
        QFETCH(double, ppm);
        std::mt19937 generator(20260929);
        for (int trial = 0; trial < 256; ++trial) {
            QVector<QPointF> points;
            const int count = 1 + generator() % 500;
            for (int i = 0; i < count; ++i) {
                const double mass = 100.0 + (generator() % 100) * 0.001;
                double intensity = static_cast<int>(generator() % 8) - 2;
                if (i % 31 == 0) intensity = -0.0;
                if (i % 37 == 0) intensity = std::numeric_limits<double>::quiet_NaN();
                if (i % 41 == 0) intensity = std::numeric_limits<double>::infinity();
                if (i % 43 == 0) intensity = -std::numeric_limits<double>::infinity();
                points.push_back({mass, intensity});
            }
            QVector<QPointF> targets = {{99.0, 1.0}, {101.0, 1.0}};
            for (int i = 0; i < 24; ++i) {
                const double mass = 100.0 + (generator() % 100) * 0.001;
                targets.push_back({mass, 1.0});
            }
            QVector<QPointF> sorted = points;
            std::sort(sorted.begin(), sorted.end(),
                      [](const QPointF &a, const QPointF &b) { return a.x() < b.x(); });
            const auto expected = originalExtraction(points, targets, ppm);
            const auto actual = MsUtils::extractPointsFromSortedPoints(sorted, targets, ppm);
            const auto compatibility = MsUtils::extractPointsFromPoints(points, targets, ppm);
            comparePoints(actual.mzFoundVsSearched, expected.mzFoundVsSearched);
            comparePoints(actual.intensityFoundVsSearched, expected.intensityFoundVsSearched);
            comparePoints(compatibility.mzFoundVsSearched, expected.mzFoundVsSearched);
            comparePoints(compatibility.intensityFoundVsSearched, expected.intensityFoundVsSearched);
            QVector<double> masses;
            for (const auto &target : targets) masses.push_back(target.x());
            for (bool removeZero : {false, true}) {
                QVector<QPointF> converted;
                for (int i = 0; i < expected.mzFoundVsSearched.size(); ++i) {
                    const double intensity = expected.intensityFoundVsSearched.at(i).x();
                    if (removeZero && (intensity < 0 || MathUtils::tZero(intensity))) continue;
                    converted.push_back({expected.mzFoundVsSearched.at(i).y(), std::max(intensity, 0.0)});
                }
                comparePoints(MsUtils::extractPointsFromSortedPoints(sorted, masses, ppm, removeZero), converted);
                comparePoints(MsUtils::extractPointsFromPoints(points, masses, ppm, removeZero), converted);
            }
        }
    }
};
QTEST_MAIN(PreparedSpectrumTests)
#include "PreparedSpectrumTests.moc"
