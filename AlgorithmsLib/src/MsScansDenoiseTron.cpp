//
// Created by anichols on 2/23/23.
//

#include "MsScansDenoiseTron.h"

Err MsScansDenoiseTron::init(const PythiaParameters &pythiaParameter) {

    ERR_INIT

    m_pythiaParameters = pythiaParameter;
    //TODO check pythia parameters are valid.

    ERR_RETURN
}

Err MsScansDenoiseTron::denoiseScansFrame(const QMap<ScanNumber, ScanPoints> &scansFrame) {

    ERR_INIT



    ERR_RETURN
}
