//
// Created by Drucifer on 11/20/2021.
//

#ifndef PYTHIACPP_PROTEINDIGESTOMATIC_H
#define PYTHIACPP_PROTEINDIGESTOMATIC_H

#include "WorkFlowsLib_Exports.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "Error.h"

#include <QStringList>


using namespace Error;

class WORKFLOWSLIB_EXPORTS IstrosLibraryBuilder{

public:

    explicit IstrosLibraryBuilder(const PythiaParameters &pythiaParameters);

    ~IstrosLibraryBuilder() = default;

    static Err exec();


private:

    PythiaParameters m_pythiaParameters;

};







#endif //PYTHIACPP_PROTEINDIGESTOMATIC_H
