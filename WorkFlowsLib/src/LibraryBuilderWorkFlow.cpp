//
// Created by anichols on 2/17/23.
//

#include "LibraryBuilderWorkFlow.h"

#include "PeptidesLibraryTron.h"

Err LibraryBuilderWorkFlow::exec(
        const PythiaParameters &pythiaParameters,
        const QString &filePath
        ) {

    ERR_INIT

    PeptidesLibraryTron peptidesLibraryTron(pythiaParameters);
    e = peptidesLibraryTron.exec(filePath, 666); ree;



    ERR_RETURN

}
