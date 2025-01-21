//
// Created by anichols on 8/24/22.
//

#ifndef PYTHIACPP_MSREADERPOINTERACC_H
#define PYTHIACPP_MSREADERPOINTERACC_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "MsReaderBase.h"

#include <QSharedPointer>


using namespace Error;

class FILEREADERSLIB_EXPORTS MsReaderPointerAcc {

public:

    MsReaderPointerAcc();

    ~MsReaderPointerAcc() = default;

    /**
    * @brief Opens a file using MsReaderPointerAcc.
    *
    * This function initializes the MsReaderPointerAcc with the specified file path (`filePath`).
    * It sets the MsReaderPointer to the corresponding MsReader instance associated with the
    * provided file path. Any encountered errors during the process are indicated by the returned Err code.
    *
    * @param filePath The file path to be opened with the MsReaderPointerAcc.
    * @return Err The error code indicating success or failure of the operation.
    *
    * @see setMsReaderPointer
    */
    Err openFile(const QString &filePath);

    /**
    * @brief Opens a file using MsReaderPointerAcc with optional column filtering.
    *
    * This function initializes the MsReaderPointerAcc based on the file type determined by
    * the file extension. If the file is of type MzML, it creates an instance of MsReaderMzML,
    * and if it is of type Parquet or Cached, it creates an instance of MsReaderParquet.
    * The MsReaderPointer is then set to the corresponding MsReader instance, and the file is
    * opened with optional column filtering based on the specified `columnToFilterBy` and
    * `filterRange`. Any encountered errors during the process are indicated by the returned Err code.
    *
    * @param filePath The file path to be opened with the MsReaderPointerAcc.
    * @param columnToFilterBy The name of the column to filter data by (optional).
    * @param filterRange The numeric range within which data is filtered (optional).
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange
    );

    /**
     * eFunctionNotImplemented
     *
     * @param filePath
     * @param columnToFilterBy
     * @return
     */
    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy
    );

    /**
    * @brief Checks if the MsReaderPointerAcc is initialized.
    *
    * This function checks if the MsReaderPointerAcc is initialized by verifying if the
    * MsReaderPointer (`ptr`) is not null. If `ptr` is not null, it further checks if
    * the underlying MsReader instance is initialized using its own `isInit` method.
    * The function returns true if both conditions are met, otherwise, it returns false.
    *
    * @return bool True if MsReaderPointerAcc is initialized, false otherwise.
    *
    */
    bool isInit();

    void setUseLazyLoading(bool useLazyLoading);

    QSharedPointer<MsReaderBase> ptr;

    bool useLazyLoading() const;


private:

    Err setMsReaderPointer(const QString &filePath);

private:

    QString m_cachedFilePath;
    QMap<MzTargetKey, FilePath> m_mzTargetKeyVsFilePathCache;
    bool m_useLazyLoading;

};


#endif //PYTHIACPP_MSREADERPOINTERACC_H
