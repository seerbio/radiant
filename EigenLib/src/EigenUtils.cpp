//
// Created by Drucifer on 11/7/2021.
//

#include "EigenUtils.h"

#include <QPointF>

Eigen::VectorX<double> EigenUtils::convertQPointFVecToEigen(
        const QVector<QPointF> &points,
        int precision,
        double valMax
        ) {

    const int vecLen = MathUtils::hashDecimal(valMax, precision);

    Eigen::VectorX<double> vec(vecLen);
    vec.setZero();

    for (const QPointF &p : points) {

        const int xHashed = MathUtils::hashDecimal(p.x(), precision);

        if (xHashed >= vecLen) {
            continue;
        }

        vec.coeffRef(xHashed) = p.y();
    }

    return vec;
}
