//
// Created by anichols on 2/18/23.
//

#ifndef PYTHIADIACPP_PARQUETREADERBASE_H
#define PYTHIADIACPP_PARQUETREADERBASE_H

#include "FileReadersLib_Exports.h"
#include "SqlUtils.h"

#include <arrow/api.h>
#include <arrow/csv/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <QDebug>

#include <iostream>


class FILEREADERSLIB_EXPORTS ParquetReaderBase {

public:

    ParquetReaderBase() = default;

    ~ParquetReaderBase() = default;

    virtual bool readFile(const std::string &fileURI) = 0;

    virtual bool writeFile(const std::string &outputFilePath) = 0;

    static bool isEqual(int a, int b);

    static std::map<std::string, int> buildSchemaMap(const std::shared_ptr<arrow::Schema> &schema);

    template<typename T, typename U>
    static bool readColumn(
            const std::shared_ptr<arrow::Table> &table,
            const std::shared_ptr<arrow::Schema> &schema,
            int columnIndex,
            std::vector<U> *output
    ) {

        if (schema->field_names().empty()) {
            return false;
        }

        std::shared_ptr<arrow::ChunkedArray> mzCol = table->column(columnIndex);
        const std::shared_ptr<arrow::Array> mzColChunks = mzCol->chunks()[0];

        auto arrowArray = std::static_pointer_cast<T>(mzColChunks);

        for (int64_t i = 0; i < mzColChunks->length(); ++i) {
            output->push_back(arrowArray->Value(i));
        }

        return true;
    }

    template<typename ArrowBuilderType, typename T>
    static arrow::Status buildParquetDataArray(
            const std::vector<T> &vec,
            std::shared_ptr<arrow::Array> *output
    ) {

        arrow::Status st;

        ArrowBuilderType arraybuilder;
        ARROW_RETURN_NOT_OK(arraybuilder.AppendValues(&(vec[0]), vec.size()));
        ARROW_ASSIGN_OR_RAISE(*output, arraybuilder.Finish());

        return st;
    }

    static arrow::Status buildParquetDataArrayFromString(
            const std::vector<std::string> &vec,
            std::shared_ptr<arrow::Array> *output
    ) {

        arrow::Status st;

        arrow::StringBuilder stringBuilder = arrow::StringBuilder();

        ARROW_RETURN_NOT_OK(stringBuilder.AppendValues(vec));
        ARROW_ASSIGN_OR_RAISE(*output, stringBuilder.Finish());

        return st;
    }

    static arrow::Status buildParquetDataArrayFromByteArray(
            const std::vector<std::string> &vec,
            std::shared_ptr<arrow::Array> *output
    ) {

        arrow::Status st;

        arrow::BinaryBuilder binaryBuilder = arrow::BinaryBuilder();

        ARROW_RETURN_NOT_OK(binaryBuilder.AppendValues(vec));
        ARROW_ASSIGN_OR_RAISE(*output, binaryBuilder.Finish());

        return st;
    }

    template<typename T>
    static std::string qVectorToBytesStdString(const QVector<T> &vec) {
        const QByteArray arr = SqlUtils::encodeBLOB<T>(vec);
        std::string arrStr = arr.toStdString();
        return arrStr;
    }

    template<typename T>
    static QVector<T> bytesStdStringToQVector(const std::string &bytesStdString) {
        const QVector<T> vec = SqlUtils::decodeBLOB<T>(QByteArray::fromStdString(bytesStdString));
        return vec;
    }

protected:

    bool checkParquetStatus();

};


#endif //PYTHIADIACPP_PARQUETREADERBASE_H
