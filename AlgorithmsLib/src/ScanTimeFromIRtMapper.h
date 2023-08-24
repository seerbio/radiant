//
// Created by anichols on 8/24/23.
//

#ifndef PYTHIADIACPP_SCANTIMEFROMIRTMAPPER_H
#define PYTHIADIACPP_SCANTIMEFROMIRTMAPPER_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <string>
#include <vector>


using namespace Error;


class ALGORITHMSLIB_EXPORTS ScanTimeFromIRtMapper {

public:

    ScanTimeFromIRtMapper();
    ~ScanTimeFromIRtMapper() = default;

    Err init(const QString &iRTRecalibrationFilePath);

    Err predictScanTime(double iRT, double *scanTime);


private:

    Err mapRT(const QVector<QPair<double, double>> &data);

private:

    QVector<double> m_coeffs;
    QVector<double> m_points;

    int m_segments;
    int m_rtSegments;
    int m_minRtPredBin;



};


#endif //PYTHIADIACPP_SCANTIMEFROMIRTMAPPER_H
