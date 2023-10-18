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

    Err init(const QVector<RTreePointData2D> &rTreeDataPoint2D) override;

    Err init(const QVector<RTreeBoxData2D> &rTreeDataPoint2D) override;

private:

    Q_DISABLE_COPY(BoostRTreePoints) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_BOOSTRTREEPOINTS_H

