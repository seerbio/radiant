//
// Created by anichols on 4/5/23.
//

#ifndef PYTHIADIACPP_PARQUETREADER_H
#define PYTHIADIACPP_PARQUETREADER_H

#include "Error.h"
#include "ErrorUtils.h"
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

    virtual Err initFromRead(const ParquetReaderInputBase &row) {return Error::eNoError;}

    void setDataMap(const QMap<QString, QVariant> &dataMap) {
        m_dataMap = dataMap;
    }

    [[nodiscard]] QMap<QString, QVariant> dataMap() const {
        return m_dataMap;
    }

    template<typename T>
    static QByteArray qVectorToQByteArray(const QVector<T> &vec) {
        return SqlUtils::encodeBLOB<T>(vec);;
    }

    template<typename T>
    static QVector<T> bytesArrayToQVector(const QByteArray &bytesString) {
        const QVector<T> vec = SqlUtils::decodeBLOB<T>(bytesString);
        return vec;
    }

    static QString joinQStringList(const QStringList &l) {
        return l.join(S_GLOBAL_SETTINGS.SEPARATOR);
    }

    template<typename T>
    static QVector<QSharedPointer<ParquetReaderInputBase>> convertInputStructToSharedPointers(
            const QVector<T> &paraquetReaderInputBaseDerivedStruct
            ) {
        QVector<QSharedPointer<ParquetReaderInputBase>> ptrs;
        for (const T &tr : paraquetReaderInputBaseDerivedStruct) {
            QSharedPointer<ParquetReaderInputBase> ptr(new T(tr));
            ptrs.push_back(ptr);
        }

        return ptrs;
    }

    template<typename T>
    static Err convertSharedPointersToInputStruct(
            const QVector<ParquetReaderInputBase> &parquetReaderInputBases,
            QVector<T> *outputStructs
    ) {

        ERR_INIT

        for (const ParquetReaderInputBase &prib : parquetReaderInputBases) {
            T strct;
            e = strct.initFromRead(prib); ree;
            outputStructs->push_back(strct);
        }

        ERR_RETURN
    }

    static bool checkIfExpectedKeysArePresent(
            const QMap<QString, QVariant> &dataMap,
            const QStringList &keysToCheck
            ) {
        const QStringList &mapKeys = dataMap.keys();
        auto keyCheckLogic = [mapKeys](const QString &s){return mapKeys.contains(s);};
        const bool allKeysPresent = std::all_of(
                keysToCheck.begin(),
                keysToCheck.end(),
                keyCheckLogic
        );

        return allKeysPresent;
    }


protected:

    QMap<QString, QVariant> m_dataMap;

};

class ParquetReader {

public:

    ParquetReader();
    ~ParquetReader();

    template<typename T>
    static Err read(
            const QString &fileURI,
            QVector<T> *readerRows
    ) {

        ERR_INIT

        readerRows->clear();

        e = ErrorUtils::fileExists(fileURI); ree;

        ParquetReader reader;

        QVector<ParquetReaderInputBase> ptrsRead;
        e = reader.readDataFromParquet(
                fileURI,
                &ptrsRead
        ); ree;

        e = ParquetReaderInputBase::convertSharedPointersToInputStruct(
                ptrsRead,
                readerRows
        ); ree;

        ERR_RETURN
    }

    template<typename T>
    static Err read(
            const QString &fileURI,
            const QString &column,
            const QPair<double, double> &range,
            QVector<T> *readerRows
            ) {

        ERR_INIT

        e = ErrorUtils::fileExists(fileURI); ree;

        readerRows->clear();

        ParquetReader reader;

        QVector<ParquetReaderInputBase> ptrsRead;
        e = reader.readDataFromParquet(
                fileURI,
                column,
                range,
                &ptrsRead
        ); ree;

        e = ParquetReaderInputBase::convertSharedPointersToInputStruct(
                ptrsRead,
                readerRows
        ); ree;

        ERR_RETURN
    }

    template<typename T>
    static Err read(
            const QString &fileURI,
            const QString &column,
            QVector<T> *readerRows
    ) {

        ERR_INIT

        readerRows->clear();

        ParquetReader reader;

        QVector<ParquetReaderInputBase> ptrsRead;
        e = reader.readDataFromParquetUniqueByColumn(
                fileURI,
                column,
                &ptrsRead
        ); ree;

        e = ParquetReaderInputBase::convertSharedPointersToInputStruct(
                ptrsRead,
                readerRows
        ); ree;

        ERR_RETURN
    }

    template<typename T>
    static Err write(
            const QVector<T> &readerRows,
            const QString &outputFilePath
    ) {

        ERR_INIT

        const QVector<QSharedPointer<ParquetReaderInputBase>> ptrs
                = ParquetReaderInputBase::convertInputStructToSharedPointers(readerRows);

        ParquetReader reader;
        e = reader.writeDataToParquet(
                outputFilePath,
                ptrs
        ); ree;

        ERR_RETURN
    }


    Err readDataFromParquetUniqueByColumn(
            const QString &parquetFilePath,
            const QString &uniqueColumn,
            QVector<ParquetReaderInputBase> *rowsRead
            );

    Err writeDataToParquet(
            const QString &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
            );

    Err readDataFromParquet(
            const QString &parquetFilePath,
            QVector<ParquetReaderInputBase> *rowsRead
            );

    Err readDataFromParquet(
            const QString &parquetFilePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange,
            QVector<ParquetReaderInputBase> *rowsRead
    );

private:

    Q_DISABLE_COPY(ParquetReader) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_PARQUETREADER_H
