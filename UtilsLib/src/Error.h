//
// Created by anichols on 11/07/2021.
//

#ifndef ERROR_H
#define ERROR_H

#include "UtilsLib_Exports.h"


#include <QDebug>
#include <QMap>

namespace Error{

    enum UTILSLIB_EXPORTS Err {
        eNoError = 0,
        eError,
        eEmptyContainerError,
        eValueError,
        eMemoryError,
        eFileError,
        eFileIncorrectTypeError,
        eNetworkError,
        eFunctionNotImplemented
    };


extern const UTILSLIB_EXPORTS QMap<Err, QString> errorMap;

#define  ERR_INIT Err e = eNoError;
#define ERR_RETURN return e;

// If error, log it and return it.
#define ree if (e) { qDebug() << QString::fromStdString(__FILE__) + "(" + QString::number(__LINE__) + ")" << errorMap.value(e); return e; }

#define rree if (e) { qDebug() << QString::fromStdString(__FILE__) + "(" + QString::number(__LINE__) + ")" << errorMap.value(e); return {e, {}}; }

#define einfo { qDebug() << QString::fromStdString(__FILE__) + "(" + QString::number(__LINE__) + ")" ;}

//if error, log it and return its object.
#define rrr(...) { qDebug() << QString::fromStdString(__FILE__) + "(" + QString::number(__LINE__) + ")" << errorMap.value(__VA_ARGS__); return __VA_ARGS__; }

//if error, return no error.
#define eee_absorb { qDebug() << QString::fromStdString(__FILE__) + "(" + QString::number(__LINE__) + ")" << errorMap.value(e); e = eNoError;}

}//NAMESPACE


#endif //ERROR_H
