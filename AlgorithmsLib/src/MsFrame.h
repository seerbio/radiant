 //
// Created by anichols on 4/8/23.
//

#ifndef PYTHIADIACPP_MSFRAME_H
#define PYTHIADIACPP_MSFRAME_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


class ALGORITHMSLIB_EXPORTS MsFrame {

public:

    friend class MsFrameTests;

    MsFrame() = default;
    ~MsFrame() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints> &scanPoints
            );

    Err preprocessMsFrame(
            bool denoise,
            bool deisotope,
            bool smooth
            );

private:

    Err denoiseFrame();
    Err deisotopeFrame();
    Err smoothFrame();

private:

    PythiaParameters m_params;
    QMap<ScanNumber, ScanPoints> m_frame;


};


#endif //PYTHIADIACPP_MSFRAME_H
