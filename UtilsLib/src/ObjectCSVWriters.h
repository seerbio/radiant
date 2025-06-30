//
// Created by andrewnichols on 8/8/24.
//

#ifndef OBJECTCSVWRITERS_H
#define OBJECTCSVWRITERS_H

#include "UtilsLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"

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

    template <typename T, typename U>
    static Err writeVectorToFile(
        const QVector<QPair<T, U>> &vec,
        const QString &filePathDestination
        ) {

            ERR_INIT

            QFile file(filePathDestination);
            e = ErrorUtils::isTrue(file.open(QIODevice::WriteOnly)); ree;

            QTextStream stream(&file);

            for (const QPair<T, U> &v : vec) {
                stream << v.first << S_GLOBAL_SETTINGS.COMMA << v.second <<  '\n';
            }

            file.close();

            ERR_RETURN
        }

    template <typename T, typename U>
    static Err readVectorFromFile(
        const QString &filePath,
        QVector<QPair<T, U>> *vec
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
            const QStringList lineSplit = line.split(S_GLOBAL_SETTINGS.COMMA);

            e = ErrorUtils::isEqual(lineSplit.size(), 2); ree;

            double output1;
            double output2;
            e = ErrorUtils::toDouble(lineSplit.front(), &output1); ree;
            e = ErrorUtils::toDouble(lineSplit.back(), &output2); ree;

            vec->push_back({static_cast<T>(output1), static_cast<U>(output2)});
        }

        ERR_RETURN
    }

    static Err writeScanPoints(
        const QMap<int, QVector<PointFF>> &scanPoints,
        const QString &filePathDestination
        );

	static Err writeScanPoints(
		const QVector<PointFF> &scanPoints,
		const QString &filePathDestination
		);

	static Err readScanPoints(
		const QString &filePathSource,
		QVector<PointFF> *scanPoints
		);

};



#endif //OBJECTCSVWRITERS_H
