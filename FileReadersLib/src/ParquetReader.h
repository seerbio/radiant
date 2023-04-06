//
// Created by anichols on 4/5/23.
//

#ifndef PYTHIADIACPP_PARQUETREADER_H
#define PYTHIADIACPP_PARQUETREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "SqlUtils.h"

#include <QDebug>

using namespace Error;

struct FILEREADERSLIB_EXPORTS ParquetReaderInputBase {

public:

    ParquetReaderInputBase() = default;
    ~ParquetReaderInputBase() = default;

    virtual QMap<QString, QVariant> map() {return {};}

    void setDataMap(const QMap<QString, QVariant> &dataMap) {
        m_dataMap = dataMap;
    }

    QMap<QString, QVariant> dataMap() const {
        return m_dataMap;
    }

    template<typename T>
    static QByteArray qVectorToQByteArray(const QVector<T> &vec) {
//        const QByteArray arr = SqlUtils::encodeBLOB<T>(vec);
//        std::string arrStr = arr.toStdString();
        return SqlUtils::encodeBLOB<T>(vec);;
    }

    template<typename T>
    static QVector<T> bytesStdStringToQVector(const std::string &bytesStdString) {
        const QVector<T> vec = SqlUtils::decodeBLOB<T>(QByteArray::fromStdString(bytesStdString));
        return vec;
    }

    static QString joinQStringList(const QStringList &l) {
        return l.join(S_GLOBAL_SETTINGS.SEPARATOR);
    }

protected:

    QMap<QString, QVariant> m_dataMap;

};

class ParquetReader {

public:

    ParquetReader();
    ~ParquetReader();

    Err writeDataToParquet(
            const QString &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
            );

    Err readDataFromParquet(
            const QString &parquetFilePath,
            QVector<ParquetReaderInputBase> *rowsRead
            );

private:

    Q_DISABLE_COPY(ParquetReader) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_PARQUETREADER_H
