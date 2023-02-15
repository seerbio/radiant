//
// Created by anichols on 11/07/2021.
//

#ifndef ERROR_H
#define ERROR_H

#include <map>
#include <iostream>

namespace Error{

    /*
    * When adding values to this enum, make sure to add the name and index
    * to std::map errorMap in the .cpp file.
    */
    enum Err {
        eNoError = 0,
        eError,
        eValueError,
        eMemoryError,
        eFileError,
        eFileIncorrectTypeError,
        eNetworkError,
        eFunctionNotImplemented,
    };


extern const std::map<Err, std::string> errorMap;

#define  ERR_INIT Err e = eNoError;
#define ERR_RETURN return e;

// If error, log it and return it.
#define ree if (e) { std::cout << __FILE__ << "LINE " << __LINE__ << " " << errorMap.at(e) << std::endl; return e; }

//if error, log it and return its object.
#define rrr(...) { std::cout << __FILE__ << "LINE " << __LINE__ << " " << errorMap.at(__VA_ARGS__) << std::endl; return __VA_ARGS__; }

//if error, return no error.
#define eee_absorb { std::cout << __FILE__ << "LINE" << __LINE__ << errorMap.at(e) << std::endl; e = eNoError;}

}//NAMESPACE


#endif //ERROR_H
