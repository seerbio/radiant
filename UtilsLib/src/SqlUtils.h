//
// Created by Drucifer on 11/21/2021.
//

#ifndef SQLUTILS_H
#define SQLUTILS_H

#include "UtilsLib_Exports.h"

#include <QByteArray>


class UTILSLIB_EXPORTS SqlUtils {

public:

    template<typename T>
    static QVector<T> decodeBLOB(const QByteArray &bytes) {

        T *ptr = (T*)bytes.constData();

        const int bytesSize = bytes.size();
        const int size = bytesSize/sizeof(T);

        QVector<T> returnArr;
        returnArr.reserve(size);

        for (int i = 0; i < size; i++) {
            returnArr.push_back(ptr[i]);
        }

        return returnArr;
    }


    template<typename T>
    static QByteArray encodeBLOB(const QVector<T> &vec) {

        QByteArray data = QByteArray::fromRawData (
                reinterpret_cast<const char*>(vec.constData()),
                sizeof(T) * vec.size()
        );

        return data;
    }


    static QByteArray fileChecksum(const QString &filePath);

};


#endif //SQLUTILS_H
