//
// Created by Drucifer on 4/27/2022.
//

#ifndef PYTHIACPP_MSUTILS_H
#define PYTHIACPP_MSUTILS_H

#include "AlgorithmsFFLib_Exports.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QFile>
#include <QPointF>

using namespace Error;

struct ALGORITHMSFFLIB_EXPORTS ExtractPoints {
    QVector<QPointF> mzFoundVsSearched;
    QVector<QPointF> intensityFoundVsSearched;
};


class ALGORITHMSFFLIB_EXPORTS MsUtils {


public:

    /**
    * @brief Extract data points from a set of points based on specified points and ppm threshold.
    *
    * This function extracts data points from a set of points based on specified target points and ppm threshold.
    *
    * @param _points The source set of points.
    * @param _pointsToExtract The target set of points to extract.
    * @param extractionPPM The PPM (parts per million) threshold for extraction.
    * @return An ExtractPoints object containing extracted points.
    */
	template<typename Point>
    static ExtractPoints extractPointsFromPoints(
            const QVector<Point> &_points,
            const QVector<QPointF> &_pointsToExtract,
            double extractionPPM
    ) {
	    ExtractPoints extractPointsOutput;
	    extractPointsOutput.mzFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});
	    extractPointsOutput.intensityFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});

	    QVector<QPointF> pointsToExtract = _pointsToExtract;
	    std::sort(pointsToExtract.begin(), pointsToExtract.end(), [](const QPointF &l, const QPointF &r){return l.x() < r.x();});

	    QVector<Point> points = _points;
	    std::sort(points.begin(), points.end(), [](const Point &l, const Point &r){return l.x() < r.x();});

	    int currentExtractionIndex = 0;
	    double extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
	    double extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

	    extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
	    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;

	    double xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
	    double xLo = extractionPointX - xPPM;
	    double xHi = extractionPointX + xPPM;

	    for (const Point &xPoint : points) {

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
	                    static_cast<double>(xPoint.y()),
	                    static_cast<double>(extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x())
	                    );
	        }
	    }

	    while (extractPointsOutput.intensityFoundVsSearched.back().y() < 0) {
	        extractPointsOutput.mzFoundVsSearched.pop_back();
	        extractPointsOutput.intensityFoundVsSearched.pop_back();
	    }

	    return extractPointsOutput;
    }

    /**
    * @brief Extracts data points from a set of points based on specified extraction points and ppm threshold.
    *
    * This function extracts data points from a set of points based on specified extraction points and ppm threshold.
    *
    * @param points The source set of points.
    * @param extractionPoints The set of extraction points for mz values.
    * @param extractionPPM The PPM (parts per million) threshold for extraction.
    * @param removeZeroPoints Flag to remove points with negative mz values.
    * @return A QVector<QPointF> containing the extracted points.
    */
    static QVector<QPointF> extractPointsFromPoints(
            const QVector<QPointF> &points,
            const QVector<double> &extractionPoints,
            double extractionPPM,
            bool removeZeroPoints = true
    );

    /**
    * @brief Writes points to a CSV file.
    *
    * This template function writes points to a CSV file. The function ensures that the destination file is not present
    * before writing the data. It uses a CSVReaderInputBase derived struct for mapping point data.
    *
    * @tparam T The type of points to be written.
    * @param points The QVector containing the points.
    * @param destFilePath The destination file path for writing the CSV.
    * @return An Err object indicating success or failure.
    */
    template<typename T>
    static Err writePointsToCSV(
            const QVector<T> &points,
            const QString &destFilePath
    ) {
        ERR_INIT

        e = ErrorUtils::isNotEmpty(points); ree

        if (QFile::exists(destFilePath)) {
            if (QFile::remove(destFilePath)) {
                qDebug() << "File deleted successfully.";
            } else {
                qDebug() << "Failed to delete file.";
            }
        } else {
            qDebug() << "File does not exist.";
        }

        struct PointRow : CSVReaderInputBase {

            double x = -1.0;
            double y = -1.0;

            QMap<QString, QVariant> map() override {
                return {
                        {"x", x},
                        {"y", y}
                };
            }
        };

        const auto transformLogic = [](const T &p){
            PointRow pr;
            pr.x = p.x();
            pr.y = p.y();
            return pr;
        };

        QVector<PointRow> pointRows;

        std::transform(
                points.begin(),
                points.end(),
                std::back_inserter(pointRows),
                transformLogic
        );

        e = CSVReader::write(pointRows, destFilePath); ree;

        ERR_RETURN
    }

};


#endif //PYTHIACPP_MSUTILS_H
