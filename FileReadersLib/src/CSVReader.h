//
// Created by anichols on 4/4/23.
//

#ifndef PYTHIADIACPP_CSVREADER_H
#define PYTHIADIACPP_CSVREADER_H

#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"


using namespace Error;


class CSVReaderInputBase {

public:
    virtual ~CSVReaderInputBase() = default;

    /**
     * @brief returns a collection of key/value pairs that store the information to be written.
     *        This map must be populated in order for a class that inherits CSVReaderInputBase to be able to write to file.
     *        Every variable of a class that inherits from CSVReaderInput base must have an entry here if it is to be written.
     *        See MsParquetReaderRow::map() for an example (Note: that MsParquetReaderRow inherits from ParquetReaderInputBase,
     *        but it is implemented the same way CSVReaderInputBase is)
     * @return a map of key/value pairs to be written.
     */
    virtual QMap<QString, QVariant> map() {return {};};

    /**
    * @brief Initializes object from a CSV row read. This method must be implemented to read from file.
    *
    * This method overrides its base version to provide functionality for initializing data from a
    * row read from a CSV data file. It verifies the presence of all expected keys in the
    * data map of the row. If all keys are present, the values corresponding to them are retrieved,
    * converted to their respective types, and stored in instance variables.
    *
    * @param row A reference to a CSVReaderInputBase object that contains the data read from a row.
    *
    * @return Returns an Err object which indicates the success or failure of the operation. If the
    * operation is successful, and all keys are present, an Err object initialized with a success
    * state is returned. If any of the keys are missing, an Err object initialized with a failure
    * state is returned.
    *
    * See MsParquetReaderRow::map() for an example (Note: that MsParquetReaderRow inherits from ParquetReaderInputBase,
    *        but it is implemented the same way CSVReaderInputBase is)
    */
    virtual Err initFromRead(const CSVReaderInputBase &row) {return Error::eNoError;}

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
    * @brief Converts a QVector of CSVReaderInputBase derived structure objects into a QVector of QSharedPointers to CSVReaderInputBase.
    *
    * This is a static template method that facilitates converting a QVector containing objects of a type derived
    * from CSVReaderInputBase into a QVector of QSharedPointers to CSVReaderInputBase.
    *
    * @param csvReaderInputBaseDerivedStruct A const reference to a QVector containing objects of a
    * type derived from CSVReaderInputBase.
    *
    * @return Returns a QVector containing QSharedPointers to CSVReaderInputBase.
    */
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

    /**
    * @brief Converts a QVector of CSVReaderInputBase types into a QVector of a derived structure type.
    *
    * This is a static template method that facilitates converting a QVector of CSVReaderInputBase types
    * into a QVector of a derived structure type by initializing each derived structure from the corresponding
    * CSVReaderInputBase type in the input QVector.
    *
    * @param csvReaderInputBases A const reference to a QVector containing CSVReaderInputBase objects.
    * @param outputStructs A pointer to a QVector where the converted objects will be stored.
    *
    * @return Returns an Err object that indicates whether the conversion operation was performed
    * successfully. An Err object initialized with a success state means the operation is successful.
    * If any error occurs during the conversion process, the function will return an Err object initialized with a failure state.
    */
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

    /**
    * @brief Converts an instance of CSVReaderInputBase type into a derived structure type.
    *
    * This is a static template method that facilitates converting an instance of CSVReaderInputBase type
    * into a derived structure type by initializing each derived structure from the corresponding
    * CSVReaderInputBase type in the input QVector.
    *
    * @param csvReaderInputBase A const reference to an instance containing CSVReaderInputBase object.
    * @param outputStruct A pointer to an instance of the converted object will be stored.
    *
    * @return Returns an Err object that indicates whether the conversion operation was performed
    * successfully. An Err object initialized with a success state means the operation is successful.
    * If any error occurs during the conversion process, the function will return an Err object initialized with a failure state.
    */
    template<typename T>
    static Err convertSharedPointersToInputStruct(
            const CSVReaderInputBase &csvReaderInputBase,
            T *outputStruct
    ) {

        ERR_INIT

        T strct;
        e = strct.initFromRead(csvReaderInputBase); ree;
        *outputStruct = strct;

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


class FILEREADERSLIB_EXPORTS CSVReader {


public:

    /**
    * @brief Writes data from a QVector of a derived structure to a CSV file.
    *
    * This static template method takes a QVector of objects derived from the CSVReaderInputBase
    * and writes the data to a file using the CSVReader class. It first converts the input QVector
    * to a QVector of shared pointers to CSVReaderInputBase types, and then uses the writeDataToCSV function
    * of the CSVReader class to write the data to the CSV file.
    *
    * @param readerRow A const reference to a QVector containing objects of a type derived from CSVReaderInputBase.
    * @param outputFilePath The path to the output file where the CSV data should be written.
    *
    * @return Returns an Err object indicating the success or failure of the operation. If the operation is successful,
    * an Err object initialized with a success state is returned, else an Err object with a failure state is returned.
    */
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

    /**
    * @brief Reads data from a CSV file into a QVector of derived structures.
    *
    * This static template method reads data from a CSV file using the CSVReader class and stores the
    * data in a QVector of objects of a type derived from CSVReaderInputBase. It starts by clearing
    * the provided QVector and then uses the readDataFromCSV function of the CSVReader class to read
    * the data. It then uses the convertSharedPointersToInputStruct function to convert the read data
    * into the required format and store it in the QVector.
    *
    * @param fileURI The URI of the CSV file to read data from.
    * @param readerRows A pointer to a QVector where the read data should be stored.
    *
    * @return Returns an Err object indicating the success or failure of the operation. If the operation
    * is successful, an Err object initialized with a success state is returned. If the file does not exist
    * or an error occurs during the reading process, an Err object with a failure state is returned.
    */
    template<typename T>
    static Err read(
            const QString &fileURI,
            const QString &deliniator,
            QVector<T> *readerRows
    ) {

        ERR_INIT

        qDebug() << "Reading filepath" << fileURI;
        e = ErrorUtils::fileExists(fileURI); ree

        readerRows->clear();

        CSVReader reader;

        QVector<CSVReaderInputBase> ptrsRead;
        e = reader.readDataFromCSV(
                fileURI,
                deliniator,
                &ptrsRead
        ); ree;

        e = CSVReaderInputBase::convertSharedPointersToInputStruct(
                ptrsRead,
                readerRows
        ); ree;

        ERR_RETURN
    }

    /**
    * @brief Writes data to a CSV file.
    *
    * This method takes a QVector of shared pointers to CSVReaderInputBase objects and writes
    * each of them to a new line in the output file. The keys from the QMap of the first
    * CSVReaderInputBase object are written as the header of the CSV file.
    *
    * @param outputFilePath The path to the output file where the CSV data should be written.
    * @param rowsToWrite A QVector of shared pointers to CSVReaderInputBase objects. Each object
    * represents a row of data to be written to the CSV file.
    *
    * @return Returns an Err object indicating the success or failure of the operation. If the
    * operation is successful, an Err object initialized with a success state is returned. If the QVector is
    * empty, an error occurs while opening the file, or the QMap of any CSVReaderInputBase object
    * has a different size from the map of the first object, an Err object with a failure state is returned.
    */
    Err writeDataToCSV(
            const QString &outputFilePath,
            const QVector<QSharedPointer<CSVReaderInputBase>> &rowsToWrite
            );

    /**
    * @brief Reads data from a CSV file.
    *
    * This method reads data from a specified CSV file line by line and stores the data
    * in a QVector of CSVReaderInputBase instances. The first line read is assumed
    * to be the header. The method then checks if the number of columns in
    * the header and each data row is equal. If the number of columns matches, each column
    * value from the row is stored in a QMap, with the corresponding header key, and
    * the map is stored in a CSVReaderInputBase instance.
    *
    * @param csvFilePath A QString representing the path to the CSV file to be read.
    * @param delineator A QString representing the character used to split the row.
    * @param readRows A pointer to a QVector of CSVReaderInputBase where the read data should be stored,
    * with each row data being encapsulated in a separate CSVReaderInputBase instance in the vector.
    *
    * @return Returns an Err object, with a success state if the operation was successful. If a file error occurs or
    * if the number of columns in the header and a minimum of one of the data rows does not match, it returns an Err object with a failure state.
    */
    Err readDataFromCSV(
            const QString &csvFilePath,
            const QString &delineator,
            QVector<CSVReaderInputBase> *readRows
            );

};


#endif //PYTHIADIACPP_CSVREADER_H
