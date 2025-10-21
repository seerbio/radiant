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
    static ExtractPoints extractPointsFromPoints(
            const QVector<QPointF> &points,
            const QVector<QPointF> &extractionPoints,
            double extractionPPM
    );

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
