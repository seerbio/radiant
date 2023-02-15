//
// Created by anichols on 11/07/2021.
//

#include "Error.h"

namespace Error{

    extern const std::map<Err, std::string> errorMap = {
            {eNoError, "NoError"},
            {eError, "Error"},
            {eValueError, "ValueError"},
            {eMemoryError, "MemoryError"},
            {eFileError, "FileError"},
            {eFileIncorrectTypeError, "FileIncorrectTypeError"},
            {eFunctionNotImplemented, "FunctionNotImplemented"},
            {eNetworkError, "NetworkError"}
    };

}//NAMESPACE
