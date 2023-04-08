//
// Created by anichols on 4/7/23.
//

#include "MsFilePreProcessorWorkFlow.h"

#include "DeisotoperTandem.h"
#include "ErrorUtils.h"
#include "MsScansDenoiseTron.h"

#include <QElapsedTimer>


Err MsFilePreProcessorWorkFlow::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}

