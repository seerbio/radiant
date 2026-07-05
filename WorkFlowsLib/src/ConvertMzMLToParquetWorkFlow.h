//
// Created by anichols on 4/7/23.
//

#ifndef PYTHIADIACPP_CONVERTMZMLTOPARQUETWORKFLOW_H
#define PYTHIADIACPP_CONVERTMZMLTOPARQUETWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "MsReaderPointerAcc.h"


using namespace Error;


class WORKFLOWSLIB_EXPORTS ConvertMzMLToParquetWorkFlow {

public:

    static Err convertMzMLToParquet(
            const QString &msDataFilePath,
            QString *outputFilePath
            );

};


#endif //PYTHIADIACPP_CONVERTMZMLTOPARQUETWORKFLOW_H
