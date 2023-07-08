//
// Created by anichols on 7/7/23.
//

#ifndef PYTHIADIACPP_MS2POINTSEXTRACTOMATIC_H
#define PYTHIADIACPP_MS2POINTSEXTRACTOMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"


class ALGORITHMSLIB_EXPORTS MS2PointsExtractomatic {

public:

    MS2PointsExtractomatic();
    ~MS2PointsExtractomatic();

    Err init(const QVector<ScanPoints> &scanPoints);

    Err findScanPoints(
            double mzMin,
            double mzMax,
            QMap<Id, ScanPoint> *foundIdVsScanPoints
            );

    Err removeScanPoints(const QMap<Id, ScanPoint> &scanPoints);

    Err updatedScanPoints(const QMap<Id, ScanPoint> &scanPoints);

private:

    Q_DISABLE_COPY(MS2PointsExtractomatic) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_MS2POINTSEXTRACTOMATIC_H
