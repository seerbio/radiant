//
// Created by anichols on 4/29/23.
//

#ifndef PYTHIADIACPP_MS1FEATUREFINDER_H
#define PYTHIADIACPP_MS1FEATUREFINDER_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"
#include "PythiaParameterReader.h"

using namespace Error;

class ALGORITHMSLIB_EXPORTS Ms1FeatureFinder {

public:

    Ms1FeatureFinder() = default;
    ~Ms1FeatureFinder() = default;

    Err init(const PythiaParameters &pythiaParameters);

    Err exec(const QString &msDataFilePath);

private:

    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MS1FEATUREFINDER_H
