//
// Created by anichols on 10/18/23.
//

#include "BoostRTreeWrapper.h"

#include "BoostRTreeAbstractFactory.h"
#include "BoostRTreePoints.h"


BoostRTreeWrapper::~BoostRTreeWrapper() {}

Err BoostRTreeWrapper::init(const QVector<RTreePointData2D> &rTreeData) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(rTreeData); ree;

    QSharedPointer<BoostRTreeAbstractFactory> rtree(new BoostRTreePoints());

    m_rTree = rtree;
    m_rTree->init(rTreeData); ree;

    ERR_RETURN
}

Err BoostRTreeWrapper::init(const QVector<RTreeBoxData2D> &rTreeData) {
    return eFunctionNotImplemented;
}

Err BoostRTreeWrapper::getPoints(
        double xMin,
        double xMax,
        double yMin,
        double yMax,
        QVector<RTreePointData2D> *vals
        ) {
    ERR_INIT
    e = m_rTree->getPoints(xMin, xMax, yMin, yMax, vals); ree;
    ERR_RETURN
}

Err BoostRTreeWrapper::getBoxes(
        double xMin,
        double xMax,
        double yMin,
        double yMax,
        QVector<RTreeBoxData2D> *vals
        ) {
    ERR_INIT
    e = m_rTree->getBoxes(xMin, xMax, yMin, yMax, vals); ree;
    ERR_RETURN
}
