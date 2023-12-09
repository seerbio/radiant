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
}
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

    Err init(const QString &iRTRecalibrationFilePath);

    Err init(const QVector<QPair<XVal, YVal>> &data);

    Err predictY(double x, double *y) const;


private:

    Err mapXtoY(const QVector<QPair<double, double>> &data);

    static Err _splineTestAcces(
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
