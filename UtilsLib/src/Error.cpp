//
// Created by anichols on 11/07/2021.
//

#include "Error.h"

namespace Error{

    extern const QMap<Err, QString> errorMap = {
            {eNoError, QStringLiteral("NoError")},
            {eError, QStringLiteral("Error")},
            {eValueError, QStringLiteral("ValueError")},
            {eMemoryError, QStringLiteral("MemoryError")},
            {eFileError, QStringLiteral("FileError")},
            {eFileIncorrectTypeError, QStringLiteral("FileIncorrectTypeError")},
            {eFunctionNotImplemented, QStringLiteral("FunctionNotImplemented")},
            {eNetworkError, QStringLiteral("NetworkError")}
    };

}//NAMESPACE
