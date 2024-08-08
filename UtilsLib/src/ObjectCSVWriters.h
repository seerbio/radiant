//
// Created by andrewnichols on 8/8/24.
//

#ifndef OBJECTCSVWRITERS_H
#define OBJECTCSVWRITERS_H

#include "UtilsLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"

#include <QFile>

using namespace Error;

class UTILSLIB_EXPORTS ObjectCSVWriters {

public:

    template <typename T>
    static Err writeVectorToFile(
        const QVector<T> &vec,
        const QString &filePathDestination
        ) {

        ERR_INIT

        QFile file(filePathDestination);
        e = ErrorUtils::isTrue(file.open(QIODevice::WriteOnly)); ree;

        QTextStream stream(&file);

        for (T v : vec) {
            stream << v << '\n';
        }

        file.close();

        ERR_RETURN
    }

    template <typename T>
    static Err readVectorFromFile(
        const QString &filePath,
        QVector<T> *vec
        ) {

        ERR_INIT

        QFile file(filePath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qDebug() << "Cannot open file for reading: " << qPrintable(file.errorString());
            rrr(eFileError);
        }

        QTextStream in(&file);
        while (!in.atEnd()) {
            QString line = in.readLine();

            double output;
            e = ErrorUtils::toDouble(line, &output); ree;

            vec->push_back(static_cast<T>(output));
        }


        ERR_RETURN
    }




};



#endif //OBJECTCSVWRITERS_H
