//
// Created by anichols on 2/18/23.
//

#ifndef PYTHIADIACPP_PARQUETREADERBASE_H
#define PYTHIADIACPP_PARQUETREADERBASE_H

#include <arrow/api.h>
#include <arrow/csv/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <QDebug>

#include <iostream>



class ParquetReaderBase {

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

protected:

    bool checkParquetStatus();

};


#endif //PYTHIADIACPP_PARQUETREADERBASE_H
