//
// Created by anichols on 4/4/23.
//

#include "CSVReader.h"

#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "QFile"


Err CSVReader::writeDataToCSV(
        const QString &outputFilePath,
        const QVector<QSharedPointer<CSVReaderBase>> &rowsToWrite
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(rowsToWrite); ree;

    QFile file(outputFilePath);
    if (file.open(QIODevice::ReadWrite)) {

        const QMap<QString, QVariant> &fisrtMap = rowsToWrite.first()->map();
        const int headerSize = fisrtMap.size();

        QTextStream stream(&file);
        stream << fisrtMap.keys().join(S_GLOBAL_SETTINGS.COMMA) + S_GLOBAL_SETTINGS.NEWLINE;

        int counter = 0;
        for (const auto &rtw : rowsToWrite) {

            if (counter++ % 10000 == 0) {
                qDebug() << "writing csv line:" << counter;
            }

            const QMap<QString, QVariant> &map = rtw->map();

            e = ErrorUtils::isEqual(map.size(), headerSize); ree;

            QStringList rtwQStringList;
            for (const QVariant &qv : map) {
                rtwQStringList.push_back(qv.toString());
            }

            stream << rtwQStringList.join(S_GLOBAL_SETTINGS.COMMA);
            stream << S_GLOBAL_SETTINGS.NEWLINE;
        }

        file.close();
        qDebug() << rowsToWrite.size() << "Written to" << outputFilePath;
        ERR_RETURN
    }

    rrr(eFileError);
}

Err CSVReader::readDataFromCSV(
        const QString &outputFilePath,
        QVector<CSVReaderBase> &readRows
        ) {
    return Error::eFunctionNotImplemented;
}
