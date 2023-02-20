//
// Created by anichols on 2/19/23.
//

#include "MsFraggerTronWorkFlow.h"

#include "MsReaderMzML.h"


MsFraggerTronWorkFlow::MsFraggerTronWorkFlow(const PythiaParameters &pythiaParameters)
: m_pythiaParameters(pythiaParameters){}

Err MsFraggerTronWorkFlow::exec(const QString &mzmLFileURI) {

    ERR_INIT

    MsReaderMzML mzMLReader;
    e = mzMLReader.openFile(mzmLFileURI); ree;




    ERR_RETURN
}
