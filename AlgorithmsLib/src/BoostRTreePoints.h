//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_BOOSTRTREEPOINTS_H
#define PYTHIADIACPP_BOOSTRTREEPOINTS_H

#include "AlgorithmsLib_Exports.h"

#include "BoostRTreeAbstractFactory.h"
#include "Error.h"
#include "GlobalSettings.h"


class ALGORITHMSLIB_EXPORTS BoostRTreePoints : public BoostRTreeAbstractFactory {

public:

    BoostRTreePoints();

    ~BoostRTreePoints() override;

    Err init(const QVector<RTreePointData2D> &rTreeDataPoint2Ds) override;

    Err init(const QVector<RTreeBoxData2D> &rTreeBoxPoint2Ds) override;

    Err getPoints(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreePointData2D> *vals
    ) override;

    Err getBoxes(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreeBoxData2D> *vals
    ) override;


private:

    Q_DISABLE_COPY(BoostRTreePoints) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_BOOSTRTREEPOINTS_H

