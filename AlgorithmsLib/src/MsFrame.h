 //
// Created by anichols on 4/8/23.
//

#ifndef PYTHIADIACPP_MSFRAME_H
#define PYTHIADIACPP_MSFRAME_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

using namespace Error;


class ALGORITHMSLIB_EXPORTS MsFrame {

public:

    friend class MsFrameTests;

    MsFrame();
    ~MsFrame() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints> &scanPoints
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints> &scanPoints,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            double collisionEnergy,
            double precursorTargetMz,
            double isoWindowLower,
            double isoWindowUpper
    );

    Err preprocessMsFrame(
            bool denoise,
            bool deisotope,
            bool smooth
            );

    [[nodiscard]] QPair<double, double> precursorMzTargetStartEnd() const;

    [[nodiscard]] UniqueMsInfoScanKey uniqueMsInfoScanKey() const;

    [[nodiscard]] int scanCount() const;

    QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints() const;

private:

    Err denoiseFrame();
    Err deisotopeFrame();
    Err smoothFrame();
    Err buildFrameIndexVsScanNumber();

private:

    PythiaParameters m_params;
    QMap<ScanNumber, ScanPoints> m_frame;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    double m_collisionEnergy;
    double m_precursorTargetMz;
    double m_isoWindowLower;
    double m_isoWindowUpper;
    QMap<FrameIndex, ScanNumber> m_frameIndexVsScanNumber;

};


#endif //PYTHIADIACPP_MSFRAME_H
