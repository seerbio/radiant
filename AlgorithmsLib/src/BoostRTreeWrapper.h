//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_BOOSTRTREEWRAPPER_H
#define PYTHIADIACPP_BOOSTRTREEWRAPPER_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

class RTreePointData2D;
class RTreeBoxData2D;

class ALGORITHMSLIB_EXPORTS BoostRTreeWrapper {

public:

    BoostRTreeWrapper() = default;
    ~BoostRTreeWrapper();

    Err init(const QVector<RTreePointData2D> &rTreeData);

    Err init(const QVector<RTreeBoxData2D> &rTreeData);


private:

    Q_DISABLE_COPY(BoostRTreeWrapper) class Private;
    const QScopedPointer<Private> d_ptr;


};


#endif //PYTHIADIACPP_BOOSTRTREEWRAPPER_H
