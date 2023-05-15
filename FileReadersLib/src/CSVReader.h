//
// Created by anichols on 4/4/23.
//

#ifndef PYTHIADIACPP_CSVREADER_H
#define PYTHIADIACPP_CSVREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"


using namespace Error;


struct CSVReaderInputBase {

public:

    virtual QMap<QString, QVariant> map() {return {};};

    virtual Err initFromRead(const CSVReaderInputBase &row) {return Error::eNoError;}

    void setDataMap(const QMap<QString, QVariant> &dataMap) {
        m_dataMap = dataMap;
    }

    [[nodiscard]] QMap<QString, QVariant> dataMap() const {
        return m_dataMap;
    }

    template<typename T>
    static QVector<QSharedPointer<CSVReaderInputBase>> convertInputStructToSharedPointers(
            const QVector<T> &csvReaderInputBaseDerivedStruct
            ) {
        QVector<QSharedPointer<CSVReaderInputBase>> ptrs;
        for (const T &tr : csvReaderInputBaseDerivedStruct) {
            QSharedPointer<CSVReaderInputBase> ptr(new T(tr));
            ptrs.push_back(ptr);
        }

        return ptrs;
    }

    template<typename T>
    static Err convertSharedPointersToInputStruct(
            const QVector<CSVReaderInputBase> &csvReaderInputBases,
            QVector<T> *outputStructs
    ) {

        ERR_INIT

        for (const CSVReaderInputBase &prib : csvReaderInputBases) {
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


class FILEREADERSLIB_EXPORTS CSVReader {


public:

    template<typename T>
    static Err write(
            const QVector<T> &readerRow,
            const QString &outputFilePath
    ) {

        ERR_INIT

        const QVector<QSharedPointer<CSVReaderInputBase>> ptrs
                = CSVReaderInputBase::convertInputStructToSharedPointers(readerRow);

        CSVReader reader;
        e = reader.writeDataToCSV(
                outputFilePath,
                ptrs
        ); ree;

        ERR_RETURN
    }

    template<typename T>
    static Err read(
            const QString &fileURI,
            QVector<T> *readerRows
    ) {

        ERR_INIT

        qDebug() << "Reading filepath" << fileURI;

        readerRows->clear();

        CSVReader reader;

        QVector<CSVReaderInputBase> ptrsRead;
        e = reader.readDataFromCSV(
                fileURI,
                &ptrsRead
        ); ree;

        e = CSVReaderInputBase::convertSharedPointersToInputStruct(
                ptrsRead,
                readerRows
        ); ree;

        ERR_RETURN
    }

    Err writeDataToCSV(
            const QString &outputFilePath,
            const QVector<QSharedPointer<CSVReaderInputBase>> &rowsToWrite
            );

    Err readDataFromCSV(
            const QString &csvFilePath,
            QVector<CSVReaderInputBase> *readRows
            );

};


#endif //PYTHIADIACPP_CSVREADER_H
