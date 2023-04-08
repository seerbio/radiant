//
// Created by anichols on 8/24/22.
//

#ifndef PYTHIACPP_MSREADERPOINTERACC_H
#define PYTHIACPP_MSREADERPOINTERACC_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"

#include <QSharedPointer>


using MsReaderPointer = QSharedPointer<MsReaderBase>;

using namespace Error;

class FILEREADERSLIB_EXPORTS MsReaderPointerFactory {

public:

    static QPair<Err, MsReaderPointer> createInstance(const QString &filePath);



};


#endif //PYTHIACPP_MSREADERPOINTERACC_H
