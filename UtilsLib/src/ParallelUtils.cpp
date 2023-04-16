//
// Created by Drucifer on 12/31/2021.
//

#include "ParallelUtils.h"

QPair<QVector<double>, QVector<double>> ParallelUtils::unZip(const QVector<QPointF> &points) {

    QVector<double> v1;
    v1.reserve(points.size());

    QVector<double> v2;
    v1.reserve(points.size());

    for (const QPointF &p : points) {
        v1.push_back(p.x());
        v2.push_back(p.y());
    }

    return {v1, v2};
}
