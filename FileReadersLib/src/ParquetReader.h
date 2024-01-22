//
// Created by anichols on 4/5/23.
//

#ifndef PYTHIADIACPP_PARQUETREADER_H
#define PYTHIADIACPP_PARQUETREADER_H

#include "FileReadersLib_Exports.h"

#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "SqlUtils.h"

#include <QFileInfo>

using namespace Error;

struct FILEREADERSLIB_EXPORTS ParquetReaderInputBase {

public:

    ParquetReaderInputBase() = default;
    ~ParquetReaderInputBase() = default;

    /**
    * @brief returns a collection of key/value pairs that store the information to be written.
    *        This map must be populated in order for a class that inherits ParquetReaderInputBase to be able to write to file.
    *        Every variable of a class that inherits from ParquetReaderInput base must have an entry here if it is to be written.
    *        See MsParquetReaderRow::map() for an example.
    * @return a map of key/value pairs to be written.
    */
    virtual QMap<QString, QVariant> map() {return {};}

    /**
    * @brief Initializes object from a Parquet row read. This method must be implemented to read from file.
    *
    * This method overrides its base version to provide functionality for initializing data from a
    * row read from a Parquet data file. It verifies the presence of all expected keys in the
    * data map of the row. If all keys are present, the values corresponding to them are retrieved,
    * converted to their respective types, and stored in instance variables.
    *
    * @param row A reference to a ParquetReaderInputBase object that contains the data read from a row.
    *
    * @return Returns an Err object which indicates the success or failure of the operation. If the
    * operation is successful, and all keys are present, an Err object initialized with a success
    * state is returned. If any of the keys are missing, an Err object initialized with a failure
    * state is returned.
    *
    * See MsParquetReaderRow::map() for an example.
    */
    virtual Err initFromRead(const ParquetReaderInputBase &row) {return Error::eNoError;}

    /**
    * @brief Sets the data map.
    *
    * This is a setter method for the m_dataMap member variable. It allows to replace the
    * current data map with a new one.
    *
    * @param dataMap A const reference to a QMap containing QString keys and QVariant values.
    */
    void setDataMap(const QMap<QString, QVariant> &dataMap) {
        m_dataMap = dataMap;
    }

    /**
    * @brief Gets the data map.
    *
    * This is a getter method for the m_dataMap member variable. It allows to retrieve
    * the current data map. As it is marked with [[nodiscard]], the caller is supposed
    * to use its result and not discard it.
    *
    * @return Returns a QMap containing QString keys and QVariant values.
    */
    [[nodiscard]] QMap<QString, QVariant> dataMap() const {
        return m_dataMap;
    }

    /**
    * @brief Converts a QVector of a specified type to a QByteArray.
    *
    * This template static function takes a QVector (`vec`) of the specified type (`T`)
    * and encodes it into a QByteArray using the `SqlUtils::encodeBLOB` function. The
    * resulting QByteArray is then returned.
    *
    * @tparam T The type of elements in the QVector.
    * @param vec The QVector to be converted.
    * @return QByteArray The converted QByteArray.
    *
    * @see QVector
    * @see QByteArray
    * @see SqlUtils::encodeBLOB
    */
    template<typename T>
    static QByteArray qVectorToQByteArray(const QVector<T> &vec) {
        return SqlUtils::encodeBLOB<T>(vec);;
    }

    /**
    * @brief Converts a QByteArray to a QVector of a specified type.
    *
    * This template static function takes a QByteArray (`bytesString`) and decodes
    * it into a QVector of the specified type (`T`) using the `SqlUtils::decodeBLOB`
    * function. The resulting QVector is then returned.
    *
    * @tparam T The type of elements in the QVector.
    * @param bytesString The QByteArray to be converted.
    * @return QVector<T> The converted QVector.
    *
    */
    template<typename T>
    static QVector<T> bytesArrayToQVector(const QByteArray &bytesString) {
        const QVector<T> vec = SqlUtils::decodeBLOB<T>(bytesString);
        return vec;
    }

    /**
    * @brief Joins a QStringList into a single QString using a specified separator.
    *
    * This static function takes a QStringList and joins its elements into
    * a single QString using the separator specified in the global settings
    * (`S_GLOBAL_SETTINGS.SEPARATOR`). The resulting QString is then returned.
    *
    * @param l The QStringList to be joined.
    * @return QString The joined QString.
    *
    */
    static QString joinQStringList(const QStringList &l) {
        return l.join(S_GLOBAL_SETTINGS.SEPARATOR);
    }

    /**
    * @brief Converts a QVector of ParquetReaderInputBase derived structure objects into a QVector of QSharedPointers to ParquetReaderInputBase.
    *
    * This is a static template method that facilitates converting a QVector containing objects of a type derived
    * from ParquetReaderInputBase into a QVector of QSharedPointers to ParquetReaderInputBase.
    *
    * @param paraquetReaderInputBaseDerivedStruct A const reference to a QVector containing objects of a
    * type derived from ParquetReaderInputBase.
    *
    * @return Returns a QVector containing QSharedPointers to ParquetReaderInputBase.
    */
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

    /**
    * @brief Converts a QVector of ParquetReaderInputBase types into a QVector of a derived structure type.
    *
    * This is a static template method that facilitates converting a QVector of ParquetReaderInputBase types
    * into a QVector of a derived structure type by initializing each derived structure from the corresponding
    * ParquetReaderInputBase type in the input QVector.
    *
    * @param ParquetReaderInputBases A const reference to a QVector containing ParquetReaderInputBase objects.
    * @param outputStructs A pointer to a QVector where the converted objects will be stored.
    *
    * @return Returns an Err object that indicates whether the conversion operation was performed
    * successfully. An Err object initialized with a success state means the operation is successful.
    * If any error occurs during the conversion process, the function will return an Err object initialized with a failure state.
    */
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

    /**
    * @brief Checks if expected keys are present in the data map.
    *
    * This is a static method that compares the keys in the data map against a list of
    * expected keys. It checks if all the expected keys are present in the data map.
    *
    * @param dataMap A const reference to a QMap containing the data map where the keys are to be checked.
    * @param keysToCheck A const reference to a QStringList containing the keys expected to be present in the data map.
    *
    * @return Returns true if all expected keys are found in the data map; otherwise, it returns false.
    */
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

    /**
    * @brief Reads data of a specified type from a Parquet file.
    *
    * This template static function reads data from a Parquet file specified by
    * the `fileURI` parameter into a QVector of the specified type (`T`). It clears
    * the provided `readerRows` vector, checks if the file exists, initializes a
    * ParquetReader, and reads data into shared pointers. The data is then converted
    * from shared pointers to input structures and stored in the `readerRows` vector.
    * Any encountered errors during the process are indicated by the returned Err code.
    *
    * @tparam T The type of data to be read from the Parquet file.
    * @param fileURI The file URI of the Parquet file.
    * @param readerRows A pointer to a QVector<T> that will be populated with the read data.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
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

    /**
    * @brief Reads data of a specified type from a Parquet file within a specified column and range.
    *
    * This template static function reads data from a Parquet file specified by the `fileURI`
    * parameter, within a specified `column` and numeric `range`, into a QVector of the specified
    * type (`T`). It checks if the file exists, initializes a ParquetReader, and reads data into
    * shared pointers based on the specified column and range. The data is then converted from shared
    * pointers to input structures and stored in the `readerRows` vector. Any encountered errors
    * during the process are indicated by the returned Err code.
    *
    * @tparam T The type of data to be read from the Parquet file.
    * @param fileURI The file URI of the Parquet file.
    * @param column The name of the column from which data is to be read.
    * @param range The numeric range within which data is to be read.
    * @param readerRows A pointer to a QVector<T> that will be populated with the read data.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
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

    /**
    * @brief Reads unique data of a specified type from a Parquet file within a specified column.
    *
    * This template static function reads unique data from a Parquet file specified by the `fileURI`
    * parameter, within a specified `column`, into a QVector of the specified type (`T`). It clears
    * the provided `readerRows` vector, initializes a ParquetReader, and reads unique data into shared
    * pointers based on the specified column. The data is then converted from shared pointers to input
    * structures and stored in the `readerRows` vector. Any encountered errors during the process are
    * indicated by the returned Err code.
    *
    * @tparam T The type of data to be read from the Parquet file.
    * @param fileURI The file URI of the Parquet file.
    * @param column The name of the column from which unique data is to be read.
    * @param readerRows A pointer to a QVector<T> that will be populated with the read unique data.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
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

    /**
    * @brief Writes data of a specified type to a Parquet file.
    *
    * This template static function writes data from a QVector of the specified type (`T`)
    * specified by the `readerRows` parameter to a Parquet file with the given `outputFilePath`.
    * It converts the input structures to shared pointers, initializes a ParquetReader, and writes
    * the data to the specified output file. Any encountered errors during the process are indicated
    * by the returned Err code.
    *
    * @tparam T The type of data to be written to the Parquet file.
    * @param readerRows The QVector containing the data to be written.
    * @param outputFilePath The file path where the Parquet file will be written.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
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

    /**
    * @brief Reads unique data from a Parquet file based on a specified column.
    *
    * This function reads unique data from a Parquet file specified by the `parquetFilePath`
    * parameter. The unique data is determined based on the specified `uniqueColumn`. The
    * results are stored in a QVector of ParquetReaderInputBase objects (`rowsRead`). The
    * function uses Arrow to read and process the Parquet file, extracting columns and converting
    * them to rows. Any encountered errors during the process are indicated by the returned Arrow::Status.
    *
    * @param parquetFilePath The file path of the Parquet file to be read.
    * @param uniqueColumn The column based on which unique data is determined.
    * @param rowsRead A pointer to a QVector<ParquetReaderInputBase> that will be populated
    *                 with the read unique data.
    * @return arrow::Status The Arrow::Status indicating success or failure of the operation.
    *
    */
    Err readDataFromParquetUniqueByColumn(
            const QString &parquetFilePath,
            const QString &uniqueColumn,
            QVector<ParquetReaderInputBase> *rowsRead
            );

    /**
    * @brief Writes data to a Parquet file.
    *
    * This function writes data to a Parquet file specified by the `outputFilePath`
    * parameter. The data to be written is provided in the `rowsToWrite` vector of
    * shared pointers to ParquetReaderInputBase. It checks for non-empty input data
    * and output file path, performs the write operation using Arrow, and verifies
    * the success of the operation. Any encountered errors during the process are
    * indicated by the returned Err code.
    *
    * @param outputFilePath The file path of the Parquet file to be written.
    * @param rowsToWrite A QVector of shared pointers to ParquetReaderInputBase
    *                    containing the data to be written.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err writeDataToParquet(
            const QString &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
            );

    /**
    * @brief Reads data from a Parquet file.
    *
    * This function reads data from a Parquet file specified by the `parquetFilePath`
    * parameter. The results are stored in a QVector of ParquetReaderInputBase objects
    * (`rowsRead`). It uses Arrow to read and process the Parquet file, extracting columns
    * and converting them to rows. Any encountered errors during the process are indicated
    * by the returned Err code.
    *
    * @param parquetFilePath The file path of the Parquet file to be read.
    * @param rowsRead A pointer to a QVector<ParquetReaderInputBase> that will be populated
    *                 with the read data.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err readDataFromParquet(
            const QString &parquetFilePath,
            QVector<ParquetReaderInputBase> *rowsRead
            );

    /**
    * @brief Reads data from a Parquet file based on a specified column and filter range.
    *
    * This function reads data from a Parquet file specified by the `parquetFilePath`
    * parameter. The data is filtered based on a specified `columnToFilterBy` and a numeric
    * `filterRange`. The results are stored in a QVector of ParquetReaderInputBase objects
    * (`rowsRead`). It uses Arrow to read and process the Parquet file, extracting columns
    * and converting them to rows. The filtered data is then appended to the `rowsRead` vector.
    * Any encountered errors during the process are indicated by the returned Arrow::Status.
    *
    * @param parquetFilePath The file path of the Parquet file to be read.
    * @param columnToFilterBy The name of the column to filter data by.
    * @param filterRange The numeric range within which data is filtered.
    * @param rowsRead A pointer to a QVector<ParquetReaderInputBase> that will be populated
    *                 with the read and filtered data.
    * @return arrow::Status The Arrow::Status indicating success or failure of the operation.
    *
    */
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
