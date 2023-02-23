#include "EigenMatrixDiskUtils.h"

#include "ErrorUtils.h"

#include <QBuffer>
#include <QDataStream>
#include <QDateTime>
#include <QFile>
#include <QVersionNumber>

#include <iostream>

namespace {

    struct FileHeader {
        quint32 magicNumber;
        QVersionNumber version;
        QDateTime dateTime;
        QString filePath;
    };

    QDataStream::Version DATASTREAM_VERSION = QDataStream::Qt_5_12;

    const QLatin1String CACHE_SUFFIX(".eigen.cache");

    /// De-serialize the Eigen::SparseMatrix to data stream
    QDataStream &operator<<(
            QDataStream &stream,
            const Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix
                            ){

        const int rows = static_cast<int>(matrix.rows());
        const int cols = static_cast<int>(matrix.cols());
        const int nonZeros = static_cast<int>(matrix.nonZeros());
        const int outerSize = static_cast<int>(matrix.outerSize());
        const int innerSize = static_cast<int>(matrix.innerSize());

        // write out the bytes??
        stream << rows << cols << nonZeros << outerSize << innerSize;

        const char *valueBytes = reinterpret_cast<const char *>(matrix.valuePtr());
        int valueBytesSize = static_cast<int>(sizeof(double) * nonZeros);
        stream.writeRawData(valueBytes, valueBytesSize);

        using StorageIndex = typename Eigen::SparseMatrix<double, Eigen::ColMajor>::StorageIndex;

        const char *outerIndexBytes = reinterpret_cast<const char *>(matrix.outerIndexPtr());
        int outerIndexBytesSize = static_cast<int>(sizeof(StorageIndex) * outerSize);
        stream.writeRawData(outerIndexBytes, outerIndexBytesSize);

        const char *innerIndexBytes = reinterpret_cast<const char *>(matrix.innerIndexPtr());
        int innerIndexBytesSize = static_cast<int>(sizeof(StorageIndex) * nonZeros);
        stream.writeRawData(innerIndexBytes, innerIndexBytesSize);

        return stream;
    }

    /// Serialize the Eigen::SparseMatrix to data stream
    /// Inspired by: https://scicomp.stackexchange.com/a/21438
    QDataStream &operator>>(
            QDataStream &stream,
            Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix
            ) {

        int rows = 0;
        int cols = 0;
        int nonZeros = 0;
        int outerSize = 0;
        int innerSize = 0;

        stream >> rows >> cols >> nonZeros >> outerSize >> innerSize;

        matrix.resize(rows, cols);
        matrix.makeCompressed();
        matrix.resizeNonZeros(nonZeros);

        char *valueBytes = reinterpret_cast<char *>(matrix.valuePtr());
        uint valueBytesSize = static_cast<uint>(sizeof(double) * nonZeros);
        stream.readRawData(valueBytes, valueBytesSize);

        using StorageIndex = typename Eigen::SparseMatrix<double, Eigen::ColMajor>::StorageIndex;

        char *outerIndexBytes = reinterpret_cast<char *>(matrix.outerIndexPtr());
        uint outerIndexBytesSize = static_cast<uint>(sizeof(StorageIndex) * outerSize);
        stream.readRawData(outerIndexBytes, outerIndexBytesSize);

        char *innerIndexBytes = reinterpret_cast<char *>(matrix.innerIndexPtr());
        uint innerIndexBytesSize = static_cast<uint>(sizeof(StorageIndex) * nonZeros);
        stream.readRawData(innerIndexBytes, innerIndexBytesSize);

        matrix.finalize();

        return stream;
    }


//    QDataStream &operator<<(QDataStream &stream, const Eigen::VectorXd &data) {

//        using VectorScalar = Eigen::VectorXd::Scalar;
//        const Eigen::Index timeIndexSize = data.size();
//        stream << timeIndexSize;

//        const char *timeIndexBytes = reinterpret_cast<const char*>(data.data());
//        int timeIndexBytesSize = static_cast<int>(timeIndexSize * sizeof(VectorScalar));
//        stream.writeRawData(timeIndexBytes, timeIndexBytesSize);

//        return stream;
//    }

//    QDataStream &operator>>(QDataStream &stream, Eigen::VectorXd &matrix)
//    {
//        using VectorScalar = Eigen::VectorXd::Scalar;
//        Eigen::Index timeIndexSize;
//        stream >> timeIndexSize;

//        matrix.resize(timeIndexSize);
//        char *timeIndexBytes = reinterpret_cast<char *>(matrix.data());
//        int timeIndexBytesSize = static_cast<int>(timeIndexSize * sizeof(VectorScalar));
//        stream.readRawData(timeIndexBytes, timeIndexBytesSize);


//        return stream;
//    }

    const QString DATETIME_FORMAT = QStringLiteral("yyyy-MM-dd HH:mm:ss");

    const quint32 HEADER_MAGIC_NUMBER = 0x70012B00; // TURBOO in https://en.wikipedia.org/wiki/Hexspeak

    const QByteArray HEADER_MARK = QByteArrayLiteral("__com.eigen.cache");
    const QByteArray VERSION_MARK = QByteArrayLiteral("eigen::Version");
    const QByteArray DATETIME_MARK = QByteArrayLiteral("eigen::DateTime");
    const QByteArray OPTIONS_MARK = QByteArrayLiteral("eigen::Options");
    const QByteArray FILE_PATH_MARK = QByteArrayLiteral("eigen::FilePath");
    const QByteArray TIME_INDEX_MARK = QByteArrayLiteral("eigen::TimeIndex");
    const QByteArray BLOB_MARK = QByteArrayLiteral("eigen::Blob");

    const QString INVALID_FILE_FORMAT = QStringLiteral("Invalid file format");

} // namespace


Err EigenMatrixDiskUtils::saveEigenMatrix(
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix,
        const QString &filePath
        ){

    ERR_INIT;

    QByteArray blob;
    QBuffer inBuffer(&blob);
    inBuffer.open(QIODevice::WriteOnly);
    QDataStream ds(&inBuffer);
    ds.setVersion(DATASTREAM_VERSION);

    // header + meta info
    ds << HEADER_MAGIC_NUMBER;
    ds << HEADER_MARK;

    ds << DATETIME_MARK;
    ds << QDateTime::currentDateTime().toString(DATETIME_FORMAT).toLatin1();

    // data block
    ds << BLOB_MARK;
    ds << matrix;

    QFile file(filePath + CACHE_SUFFIX);
    const bool isFileOpened = file.open(QIODevice::WriteOnly);
    e = ErrorUtils::isTrue(isFileOpened); ree;


    long long size = file.write(blob);
    e = ErrorUtils::isNotEqual(size, static_cast<long long>(-1)); ree;

    inBuffer.close();
    file.close();


    ERR_RETURN;
}


Err EigenMatrixDiskUtils::loadEigenMatrix(const QString &filePath, Eigen::SparseMatrix<double, Eigen::ColMajor> *matrix) {
    ERR_INIT;

    QFile file(filePath + CACHE_SUFFIX);
    file.open(QIODevice::ReadOnly);

    QByteArray data = file.readAll();
    file.close();

    QBuffer inBuffer(&data);
    inBuffer.open(QIODevice::ReadOnly);
    QDataStream ds(&inBuffer);
    ds.setVersion(DATASTREAM_VERSION);

    FileHeader header;

    ds >> header.magicNumber;
    e = ErrorUtils::isEqual(static_cast<double>(header.magicNumber),
                            static_cast<double>(HEADER_MAGIC_NUMBER)); ree;

    QByteArray headerMark;
    ds >> headerMark;
    if (headerMark != HEADER_MARK) {
        rrr(eError);
    }

    QByteArray dateTimeMark;
    ds >> dateTimeMark;
    if (dateTimeMark != DATETIME_MARK) {
        rrr(eError);
    }

    QByteArray datetime;
    ds >> datetime;
    QString dateTimeStr = QString::fromLatin1(datetime);
    header.dateTime = QDateTime::fromString(dateTimeStr, DATETIME_FORMAT);


    QByteArray blobMark;
    ds >> blobMark;

    ds >> *matrix;

    ERR_RETURN;
}


Err EigenMatrixDiskUtils::saveEigenMatrix(const Eigen::MatrixXd &matrix, const QString &filePath) {
    ERR_INIT;

    const Eigen::SparseMatrix<double, Eigen::ColMajor> m = matrix.sparseView();
    e = EigenMatrixDiskUtils::saveEigenMatrix(m, filePath); ree;

    ERR_RETURN;
}


Err EigenMatrixDiskUtils::loadEigenMatrix(const QString &filePath, Eigen::MatrixXd *matrix) {
    ERR_INIT;

    Eigen::SparseMatrix<double, Eigen::ColMajor> m;
    e = EigenMatrixDiskUtils::loadEigenMatrix(filePath, &m); ree;

    *matrix = m;

    ERR_RETURN;
}


