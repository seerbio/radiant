//
// Created by anichols on 8/24/23.
//

#ifndef PYTHIADIACPP_XYMAPPERMATIC_H
#define PYTHIADIACPP_XYMAPPERMATIC_H

#include "AlgorithmsFFLib_Exports.h"

#include "CSVReader.h"
#include "GlobalSettings.h"
#include "Error.h"

#include <string>
#include <vector>


using namespace Error;

namespace XYMappermaticRowNamespace {

    const QString X = QStringLiteral("x");
    const QString Y = QStringLiteral("y");

    const QStringList keysToCheck = {
            X,
            Y
    };
}//namespace

/**
 * See CSVReaderInputBase for documentation
 */
struct XYMappermaticRow: public CSVReaderInputBase {

    double x = -1.0;
    double y = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace XYMappermaticRowNamespace;

        return {
                {X, x},
                {Y, y}
        };
    }

    Err initFromRead(const CSVReaderInputBase &row) override {

        using namespace XYMappermaticRowNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = CSVReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        if (!allKeysPresent) {
            qDebug() << dataMap;
            qDebug() << keysToCheck;
        }

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        x = dataMap.value(X).toDouble();
        y = dataMap.value(Y).toDouble();

        ERR_RETURN
    }
};


class ALGORITHMSFFLIB_EXPORTS XYMappermatic {

public:

    friend class XYMappermaticTests;

    XYMappermatic();
    ~XYMappermatic() = default;

    /**
    * @brief Initializes the XYMappermatic with iRT recalibration data.
    *
    * This method initializes the XYMappermatic by reading iRT recalibration data from the specified file.
    *
    * @param iRTRecalibrationFilePath The file path of the iRT recalibration data.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err init(const QString &iRTRecalibrationFilePath);

    /**
    * @brief Initializes the XYMappermatic with data.
    *
    * This method initializes the XYMappermatic with the given data, sorting it by X values and mapping X to Y values.
    *
    * @param _data The data containing X and Y values.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err init(const QVector<QPair<XVal, YVal>> &data);

    Err setBinning(int val);

    /**
    * @brief Predicts Y value for a given X using the initialized spline.
    *
    * This method predicts the Y value for a given X using the initialized spline coefficients and control points.
    *
    * @param x The X value for prediction.
    * @param y Pointer to store the predicted Y value.
    * @return An error code indicating the success or failure of the prediction.
    */
    Err predictY(double x, double *y) const;


private:

    Err mapXtoY(const QVector<QPair<double, double>> &data);

    static Err _splineTestAccess(
            const QVector<QPair<double, double>> &data,
            int segments,
            QVector<double> *coeffs,
            QVector<double> *points
            );

private:

    QVector<double> m_coeffs;
    QVector<double> m_points;

    int m_segments;
    int m_xSegments;
    int m_minXPredBin;

};


#endif //PYTHIADIACPP_XYMAPPERMATIC_H
