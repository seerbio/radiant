//
// Created by anichols on 2/23/23.
//

#ifndef PYTHIADIACPP_MSSCANSDENOISETRON_H
#define PYTHIADIACPP_MSSCANSDENOISETRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "FeatureFinderHill.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

using namespace Error;


class ALGORITHMSLIB_EXPORTS MsScansDenoiseTron {

public:

    MsScansDenoiseTron() = default;
    ~MsScansDenoiseTron() = default;

    Err init(const PythiaParameters &pythiaParameter);

    Err denoiseScansFrame(const QMap<ScanNumber, ScanPoints> &scansFrame,
                          QMap<ScanNumber, ScanPoints> *scansFrameDenoised
                          );

private:

    PythiaParameters m_pythiaParameters;
    FeatureFinderParameters m_ffpParams;
    FeatureFinderHillBuilder m_featureFinderHillBuilder;

};


#endif //PYTHIADIACPP_MSSCANSDENOISETRON_H
