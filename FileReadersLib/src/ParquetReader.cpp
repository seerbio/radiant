//
// Created by anichols on 4/5/23.
//

#include <iostream>

#include <arrow/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

//NOTE that these have to go after parquet arrow import or there will be a compile error.
#include "ParquetReader.h"


class Q_DECL_HIDDEN ParquetReader::Private
{
public:

    Private();
    ~Private();

    Err writeDataToParquet(
            const QString &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
    );

    //TODO add streaming version of write parquet
    // multiple tables iteratively, see parquet::arrow::FileWriter. In WriteTable() write.h parquet

    Err readDataFromParquet(
            const QString &parquetFilePath,
            QVector<ParquetReaderInputBase> *rowsRead
    );

    arrow::Status readDataFromParquet(
            const QString &parquetFilePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange,
            QVector<ParquetReaderInputBase> *rowsRead
    );

    static arrow::Status readDataFromParquetUniqueByColumn(
            const QString &parquetFilePath,
            const QString &uniqueColumn,
            QVector<ParquetReaderInputBase> *rowsRead
    );

};

ParquetReader::Private::Private() {}

ParquetReader::Private::~Private() {}

namespace {

    QMap<QString, QVector<QVariant>> splitMapsToColumns(
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
            ) {

        QMap<QString, QVector<QVariant>> mapsAsColumns;

        for (const QSharedPointer<ParquetReaderInputBase> & prib : rowsToWrite) {

            const QMap<QString, QVariant> &map = prib->map();

            for (auto it = map.begin(); it != map.end(); it++) {

                const QString &columnName = it.key();
                const QVariant &value = it.value();

                mapsAsColumns[columnName].push_back(value);
            }

        }

        return mapsAsColumns;
    }

    template<typename ArrowBuilderType, typename T>
    arrow::Status buildParquetDataArray(
            const std::vector<T> &vec,
            std::shared_ptr<arrow::Array> *output
    ) {

        arrow::Status st;

        ArrowBuilderType arraybuilder;
        ARROW_RETURN_NOT_OK(arraybuilder.AppendValues(&(vec[0]), vec.size()));
        ARROW_ASSIGN_OR_RAISE(*output, arraybuilder.Finish());

        return st;
    }

    template<typename ArrowBuilderType>
    arrow::Status buildParquetDataArrayFromString(
            const std::vector<std::string> &vec,
            std::shared_ptr<arrow::Array> *output
    ) {
        arrow::Status st;

        ArrowBuilderType stringOrBinaryArraybuilder;
        ARROW_RETURN_NOT_OK(stringOrBinaryArraybuilder.AppendValues(vec));
        ARROW_ASSIGN_OR_RAISE(*output, stringOrBinaryArraybuilder.Finish());

        return st;
    }

    arrow::Status convertQVarianVecToArrowArray(
            const QVector<QVariant> &vec,
            std::shared_ptr<arrow::Array> *arrowArr
                    ) {

        arrow::Status st;

        const QVariant &vecFirst = vec.first();

        if (
            vecFirst.typeName() == QStringLiteral("int") ||
            vecFirst.typeName() == QStringLiteral("bool")
            ) {

            std::vector<std::int64_t> stdVec;
            for (const QVariant &var : vec) {

                if (vecFirst.typeName() == QStringLiteral("bool")) {
                    stdVec.push_back(static_cast<int>(var.toBool()));
                    continue;
                }

                stdVec.push_back(var.toInt());
            }

            st = buildParquetDataArray<arrow::Int64Builder>(stdVec, arrowArr);
            if (!st.ok()){
                einfo
                return st;
            }
        }

        else if (vecFirst.typeName() == QStringLiteral("float")) {

            std::vector<float> stdVec;
            for (const QVariant &var : vec) {
                stdVec.push_back(var.toFloat());
            }

            st = buildParquetDataArray<arrow::FloatBuilder>(stdVec, arrowArr);
            if (!st.ok()){
                einfo
                return st;
            }
        }

        else if (vecFirst.typeName() == QStringLiteral("double")) {

            std::vector<double> stdVec;
            for (const QVariant &var : vec) {
                stdVec.push_back(var.toDouble());
            }

            st = buildParquetDataArray<arrow::DoubleBuilder>(stdVec, arrowArr);
            if (!st.ok()){
                einfo
                return st;
            }
        }

        else if (
                vecFirst.typeName() == QStringLiteral("QChar") ||
                vecFirst.typeName() == QStringLiteral("QString") ||
                vecFirst.typeName() == QStringLiteral("QByteArray")
                ) {

            std::vector<std::string> stdVec;
            for (const QVariant &var : vec) {

                if (vecFirst.typeName() == QStringLiteral("QByteArray")) {
                    const std::string arrStr = var.toByteArray().toStdString();
                    stdVec.push_back(arrStr);
                    continue;
                }

                const std::string arrStr = var.toString().toStdString();
                stdVec.push_back(arrStr);
            }

            if (vecFirst.typeName() == QStringLiteral("QByteArray")) {

                st = buildParquetDataArrayFromString<arrow::BinaryBuilder>(stdVec, arrowArr);
            }
            else {
                st = buildParquetDataArrayFromString<arrow::StringBuilder>(stdVec, arrowArr);
            }
            if (!st.ok()){
                einfo
                return st;
            }
        }

        if (!st.ok()) {
            einfo
            return st;
        }

        return st;
    }

    arrow::Status buildParquetColumns(
            const QMap<QString, QVector<QVariant>> &mapColumns,
            std::vector<std::shared_ptr<arrow::Array>> *arrowColumns
            ) {

        arrow::Status st;

        std::vector<std::shared_ptr<arrow::Array>> columns;

        for (auto it = mapColumns.begin(); it != mapColumns.end(); it++) {

            const QVector<QVariant> &vec = it.value();

            std::shared_ptr<arrow::Array> arrowArr;
            st = convertQVarianVecToArrowArray(
                    vec,
                    &arrowArr
                    );

            if (!st.ok()) {
                einfo
                return st;
            }

            arrowColumns->push_back(arrowArr);
        }

        return st;
    }

    std::shared_ptr<arrow::Schema> buildArrowSchema(const QSharedPointer<ParquetReaderInputBase> &row) {

        const QMap<QString, QVariant> &rowMapExample = row->map();


        std::vector<std::shared_ptr<arrow::Field>> fields;
        for (auto it = rowMapExample.begin(); it != rowMapExample.end(); it++) {

            const QString &colName = it.key();
            const QVariant &var = it.value();
            const QString &typeName = it->typeName();

            if (typeName == "double") {
                std::shared_ptr<arrow::Field> field = arrow::field(colName.toStdString(), arrow::float64());
                fields.push_back(field);
            }

            else if (typeName == "float") {
                std::shared_ptr<arrow::Field> field = arrow::field(colName.toStdString(), arrow::float32());
                fields.push_back(field);
            }

            else if (typeName == "int" || typeName == "bool") {
                std::shared_ptr<arrow::Field> field = arrow::field(colName.toStdString(), arrow::int64());
                fields.push_back(field);
            }

            else if (
                    typeName == "QChar" ||
                    typeName == "QString" ||
                    typeName == "QByteArray"
                    ) {

                if (typeName == "QByteArray") {
                    std::shared_ptr<arrow::Field> field = arrow::field(colName.toStdString(), arrow::binary());
                    fields.push_back(field);
                }

                else {
                    std::shared_ptr<arrow::Field> field = arrow::field(colName.toStdString(), arrow::utf8());
                    fields.push_back(field);
                }
            }
        }

        std::shared_ptr<arrow::Schema> schema = arrow::schema(fields);
        return schema;
    }

    arrow::Status writeFileLogic(
            const std::string &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
    ) {

        arrow::Status st;

        const QMap<QString, QVector<QVariant>> mapColumns  = splitMapsToColumns(rowsToWrite);

        std::vector<std::shared_ptr<arrow::Array>> parquetColumns;
        st = buildParquetColumns(
                mapColumns,
                &parquetColumns
                );

        if (!st.ok()) {
            einfo
            return st;
        }

        std::shared_ptr<arrow::Schema> schema = buildArrowSchema(rowsToWrite.first());

        std::shared_ptr<arrow::Table> table;
        table = arrow::Table::Make(schema, parquetColumns);

        std::shared_ptr<arrow::io::FileOutputStream> outfile;

        ARROW_ASSIGN_OR_RAISE(
                outfile,
                arrow::io::FileOutputStream::Open(outputFilePath)
        );
        if (!st.ok()){
            return st;
        }

        PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
                *table,
                arrow::default_memory_pool(),
                outfile
        ));

        return st;
    }

}//namespace
Err ParquetReader::Private::writeDataToParquet(
        const QString &outputFilePath,
        const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(rowsToWrite); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;

    arrow::Status st = writeFileLogic(
            outputFilePath.toStdString(),
            rowsToWrite
    );

    e = ErrorUtils::isTrue(st.ok()); ree;
    qDebug() << "File written to:" << outputFilePath;

    ERR_RETURN
}

namespace {

    std::map<std::string, int> buildSchemaMap(const std::shared_ptr<arrow::Schema> &schema) {

        std::string schemaColumns = schema->ToString();
        const std::vector<std::string> schemaColumnsSplit = schema->field_names();

        std::map<std::string, int> schemaMap;

//        qDebug() << "PARQUET SCHEMA";
//        qDebug() << QString::fromStdString(schemaColumns);
        for (int i = 0; i < schemaColumnsSplit.size(); i++) {
//            qDebug() << QString::fromStdString(schemaColumnsSplit.at(i));
            schemaMap.emplace(schemaColumnsSplit.at(i), i);
        }

        return schemaMap;
    }

    bool readColumnLogic(
            const std::shared_ptr<arrow::Array> &colChunks,
            QVector<QVariant> *outputVec
            ) {

        const QString typeName = QString::fromStdString(colChunks->type()->ToString());

        if (
            typeName == QStringLiteral("int64") ||
            typeName == QStringLiteral("bool")
            ) {

            auto arrowArray = std::static_pointer_cast<arrow::Int64Array>(colChunks);
            for (int64_t i = 0; i < colChunks->length(); ++i) {
                outputVec->push_back(static_cast<int>(arrowArray->Value(i)));
            }
        }

        else if (typeName == QStringLiteral("float")) {

            auto arrowArray = std::static_pointer_cast<arrow::FloatArray>(colChunks);
            for (int64_t i = 0; i < colChunks->length(); ++i) {
                outputVec->push_back(static_cast<float >(arrowArray->Value(i)));
            }
        }

        else if (typeName == QStringLiteral("double")) {

            auto arrowArray = std::static_pointer_cast<arrow::DoubleArray>(colChunks);
            for (int64_t i = 0; i < colChunks->length(); ++i) {
                outputVec->push_back(static_cast<double>(arrowArray->Value(i)));
            }
        }

        else if (
                typeName == QStringLiteral("string") ||
                typeName == QStringLiteral("binary")
                ) {

            if (typeName == QStringLiteral("binary")) {
                auto arrowArray = std::static_pointer_cast<arrow::BinaryArray>(colChunks);
                for (int64_t i = 0; i < colChunks->length(); ++i) {
                    outputVec->push_back(QByteArray::fromStdString(
                            static_cast<std::string>(arrowArray->Value(i)))
                    );
                }
            }

            else {
                auto arrowArray = std::static_pointer_cast<arrow::StringArray>(colChunks);
                for (int64_t i = 0; i < colChunks->length(); ++i) {

                    const QString convertedVal
                            = QString::fromStdString(static_cast<std::string>(arrowArray->Value(i)));

                    if (convertedVal.size() == 1) {

                        const QChar qChar = convertedVal[0];
                        outputVec->push_back(qChar);
                        continue;
                    }

                    outputVec->push_back(convertedVal);
                }
            }
        }

        return true;
    }

    bool readColumn(
            const std::shared_ptr<arrow::RecordBatch> &recordBatch,
            const std::shared_ptr<arrow::Schema> &schema,
            int columnIndex,
            QVector<QVariant> *outputVec
            ) {

        if (schema->field_names().empty()) {
            return false;
        }

        const std::shared_ptr<arrow::Array> colChunks = recordBatch->column(columnIndex);
        return readColumnLogic(colChunks, outputVec);
    }

    bool readColumn(
            const std::shared_ptr<arrow::Table> &table,
            const std::shared_ptr<arrow::Schema> &schema,
            int columnIndex,
            QVector<QVariant> *outputVec
    ) {

        if (schema->field_names().empty()) {
            return false;
        }

        std::shared_ptr<arrow::ChunkedArray> col = table->column(columnIndex);
        const std::shared_ptr<arrow::Array> colChunks = col->chunks()[0];
        return readColumnLogic(colChunks, outputVec);
    }

    Err extractColumns(
            const std::shared_ptr<arrow::Table> &table,
            QMap<QString, QVector<QVariant>> *columnsMap
            ) {

        ERR_INIT

        const std::shared_ptr<arrow::Schema> schema = table->schema();
        const std::map<std::string, int> schemaMap = buildSchemaMap(schema);

        bool columnChecker = true;
        for (auto it = schemaMap.begin(); it != schemaMap.end(); it++) {

                const QString colName = QString::fromStdString(it->first);
                const int colIndex = it->second;

                QVector<QVariant> columnValsVec;
                columnChecker = readColumn(
                        table,
                        schema,
                        colIndex,
                        &columnValsVec
                );
                if (!columnChecker) {
                    rrr(eError);
                }

                columnsMap->insert(colName, columnValsVec);
        }

        ERR_RETURN
    }

    Err columnsMapToRowsMap(
            const QMap<QString, QVector<QVariant>> &columnsMap,
            QVector<ParquetReaderInputBase> *rowsRead
            ) {

        ERR_INIT

        const int colSize = columnsMap.first().size();

        for (const QVector<QVariant> &v : columnsMap) {
            e = ErrorUtils::isEqual(v.size(), colSize);
        }

        for (int i = 0; i < colSize; i++) {

            QMap<QString, QVariant> rowMap;

            for (auto it = columnsMap.begin(); it != columnsMap.end(); it++) {

                const QString &colName = it.key();
                const QVector<QVariant> &var = it.value();

                rowMap.insert(colName, var.at(i));
            }

            ParquetReaderInputBase prib;
            prib.setDataMap(rowMap);

            rowsRead->push_back(prib);
        }

        ERR_RETURN
    }

}//namespace
Err ParquetReader::Private::readDataFromParquet(
        const QString &parquetFilePath,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    arrow::Status st;

    arrow::MemoryPool* pool = arrow::default_memory_pool();
    arrow::fs::LocalFileSystem fileSystem;

    std::shared_ptr<arrow::io::RandomAccessFile> input
            = fileSystem.OpenInputFile(parquetFilePath.toStdString()).ValueOrDie();

    std::unique_ptr<parquet::arrow::FileReader> arrowReader;
    st = parquet::arrow::OpenFile(input, pool, &arrowReader);
    if (!st.ok()) {
        return Error::eFileError;
    }

    std::shared_ptr<arrow::Table> table;
    st = arrowReader->ReadTable(&table);
    if (!st.ok()) {
        return Error::eError;
    }

    QMap<QString, QVector<QVariant>> columnsMap;
    e = extractColumns(
            table,
            &columnsMap
            ); ree;

    e = columnsMapToRowsMap(
            columnsMap,
            rowsRead
            ); ree;

    ERR_RETURN
}

namespace {

    Err extractColumns(
            const std::shared_ptr<arrow::RecordBatch> &recordBatch,
            QMap<QString, QVector<QVariant>> *columnsMap
    ) {

        ERR_INIT

        const std::shared_ptr<arrow::Schema> schema = recordBatch->schema();
        const std::map<std::string, int> schemaMap = buildSchemaMap(schema);

        bool columnChecker = true;
        for (auto it = schemaMap.begin(); it != schemaMap.end(); it++) {

            const QString colName = QString::fromStdString(it->first);
            const int colIndex = it->second;

            QVector<QVariant> columnValsVec;
            columnChecker = readColumn(
                    recordBatch,
                    schema,
                    colIndex,
                    &columnValsVec
            );
            if (!columnChecker) {
                rrr(eError);
            }

            columnsMap->insert(colName, columnValsVec);
        }

        ERR_RETURN
    }

    Err columnsMapToRowsMap(
            const QMap<QString, QVector<QVariant>> &columnsMap,
            const QString &columnToFilterBy,
            const QPair<double, double> &range,
            QVector<ParquetReaderInputBase> *rowsRead
            ) {

        ERR_INIT

        e = ErrorUtils::isBelowThreshold(
                range.first,
                range.second,
                ErrorUtilsParam::IncludeThreshold
                ); ree;

        const int colSize = columnsMap.first().size();

        for (const QVector<QVariant> &v : columnsMap) {
            e = ErrorUtils::isEqual(v.size(), colSize);
        }

        for (int i = 0; i < colSize; i++) {

            QMap<QString, QVariant> rowMap;

            for (auto it = columnsMap.begin(); it != columnsMap.end(); it++) {

                const QString &colName = it.key();
                const QVector<QVariant> &var = it.value();

                rowMap.insert(colName, var.at(i));
            }

            ParquetReaderInputBase prib;
            prib.setDataMap(rowMap);

            e = ErrorUtils::isTrue(rowMap.contains(columnToFilterBy)); ree;
            const double &var = rowMap.value(columnToFilterBy).toDouble();

            if (!(range.first <= var && var <= range.second)) {
                continue;
            }

            rowsRead->push_back(prib);
        }

        ERR_RETURN
    }


    arrow::Status initArrowReaderBuilder(
            const QString &parquetFilePath,
            parquet::arrow::FileReaderBuilder *readerBuilder
            ) {

        arrow::Status st;

        arrow::MemoryPool* pool = arrow::default_memory_pool();
        arrow::fs::LocalFileSystem fileSystem;

        auto readerProperties = parquet::ReaderProperties(pool);
        readerProperties.set_buffer_size(4096 * 4);
        readerProperties.enable_buffered_stream();

        auto arrowReaderProps = parquet::ArrowReaderProperties();
        arrowReaderProps.set_batch_size(128 * 1024);  // default 64 * 1024

        st = readerBuilder->OpenFile(
                parquetFilePath.toStdString(),
                /*memory_map=*/false,
                readerProperties
        );
        if (!st.ok()) {
            return st;
        }

        readerBuilder->memory_pool(pool);
        readerBuilder->properties(arrowReaderProps);

        return st;
    }

}//namespace
arrow::Status ParquetReader::Private::readDataFromParquet(
        const QString &parquetFilePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    arrow::Status st;

    parquet::arrow::FileReaderBuilder readerBuilder;
    st = initArrowReaderBuilder(
            parquetFilePath,
            &readerBuilder
            );
    if (!st.ok()) {
        return st;
    }

    std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
    ARROW_ASSIGN_OR_RAISE(arrow_reader, readerBuilder.Build());

    std::shared_ptr<::arrow::RecordBatchReader> rbReader;
    st = arrow_reader->GetRecordBatchReader(&rbReader);
    if (!st.ok()) {
        return st;
    }

    for (arrow::Result<std::shared_ptr<arrow::RecordBatch>> maybeBatch : *rbReader) {

        std::shared_ptr<arrow::RecordBatch> batch = maybeBatch.ValueOrDie();

        QMap<QString, QVector<QVariant>> columnsMap;
        extractColumns(batch, &columnsMap);

        QVector<ParquetReaderInputBase> rowsReadBatch;

        e = columnsMapToRowsMap(
                columnsMap,
                columnToFilterBy,
                filterRange,
                &rowsReadBatch
        );

        if (e != eNoError) {
            st = arrow::Status::Invalid("Invalid input");
            return st;
        }

        rowsRead->append(rowsReadBatch);
    }

    return st;
}

namespace {

    Err columnsMapToRowsMap(
            const QMap<QString, QVector<QVariant>> &columnsMap,
            const QString &columnToFilterBy,
            QVector<ParquetReaderInputBase> *rowsRead
    ) {

        ERR_INIT

        const int colSize = columnsMap.first().size();

        for (const QVector<QVariant> &v : columnsMap) {
            e = ErrorUtils::isEqual(v.size(), colSize);
        }

        QMap<QVariant, bool> uniques;

        for (int i = 0; i < colSize; i++) {

            QMap<QString, QVariant> rowMap;

            for (auto it = columnsMap.begin(); it != columnsMap.end(); it++) {

                const QString &colName = it.key();
                const QVector<QVariant> &var = it.value();

                rowMap.insert(colName, var.at(i));
            }

            ParquetReaderInputBase prib;
            prib.setDataMap(rowMap);

            e = ErrorUtils::isTrue(rowMap.contains(columnToFilterBy)); ree;
            const QVariant &var = rowMap.value(columnToFilterBy);

            if (uniques.value(var)) {
                continue;
            }

            uniques.insert(var, true);
            rowsRead->push_back(prib);
        }

        ERR_RETURN
    }
}//namespace
arrow::Status ParquetReader::Private::readDataFromParquetUniqueByColumn(
        const QString &parquetFilePath,
        const QString &uniqueColumn,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    arrow::Status st;

    parquet::arrow::FileReaderBuilder readerBuilder;
    st = initArrowReaderBuilder(
            parquetFilePath,
            &readerBuilder
    );
    if (!st.ok()) {
        return st;
    }

    std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
    ARROW_ASSIGN_OR_RAISE(arrow_reader, readerBuilder.Build());

    std::shared_ptr<::arrow::RecordBatchReader> rbReader;
    st = arrow_reader->GetRecordBatchReader(&rbReader);
    if (!st.ok()) {
        return st;
    }

    for (arrow::Result<std::shared_ptr<arrow::RecordBatch>> maybeBatch : *rbReader) {

        std::shared_ptr<arrow::RecordBatch> batch = maybeBatch.ValueOrDie();

        QMap<QString, QVector<QVariant>> columnsMap;
        extractColumns(batch, &columnsMap);

        QVector<ParquetReaderInputBase> rowsReadBatch;

        e = columnsMapToRowsMap(
                columnsMap,
                uniqueColumn,
                &rowsReadBatch
        );

        if (e != eNoError) {
            st = arrow::Status::Invalid("Invalid input");
            return st;
        }

        rowsRead->append(rowsReadBatch);
    }

    return st;
}

/////////////////////////////////////////////////////
// END PRIVATE //////////////////////////////////////
/////////////////////////////////////////////////////

ParquetReader::ParquetReader() : d_ptr(new Private()) {}

ParquetReader::~ParquetReader() {}

Err ParquetReader::writeDataToParquet(
        const QString &outputFilePath,
        const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
        ) {

    ERR_INIT

    e = d_ptr->writeDataToParquet(
            outputFilePath,
            rowsToWrite
            ); ree;

    qDebug() << "PQR file saved to" << outputFilePath;

    ERR_RETURN
}

Err ParquetReader::readDataFromParquet(
        const QString &parquetFilePath,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    e = d_ptr->readDataFromParquet(
            parquetFilePath,
            rowsRead
            ); ree;

    ERR_RETURN
}

Err ParquetReader::readDataFromParquet(
        const QString &parquetFilePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    arrow::Status st = d_ptr->readDataFromParquet(
            parquetFilePath,
            columnToFilterBy,
            filterRange,
            rowsRead
    );

    e = ErrorUtils::isTrue(st.ok()); ree;

    ERR_RETURN
}

Err ParquetReader::readDataFromParquetUniqueByColumn(
        const QString &parquetFilePath,
        const QString &uniqueColumn,
        QVector<ParquetReaderInputBase> *rowsRead
        ) {

    ERR_INIT

    arrow::Status st;

    st = d_ptr->readDataFromParquetUniqueByColumn(
            parquetFilePath,
            uniqueColumn,
            rowsRead
    );

    e = ErrorUtils::isTrue(st.ok()); ree;


    ERR_RETURN

}
