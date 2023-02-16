//
// Created by anichols on 2/6/23.
//

#include "MsParquetReader.h"

//#include "ErrorUtils.h"
//#include "StringUtils.h"

#include <arrow/api.h>
#include <arrow/csv/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <iostream>


namespace {

    arrow::Status GenInitialFile() {

        // Make a couple 8-bit integer arrays and a 16-bit integer array -- just like
        // basic Arrow example.
        arrow::Int8Builder int8builder;
        int8_t days_raw[5] = {1, 12, 17, 23, 28};
        ARROW_RETURN_NOT_OK(int8builder.AppendValues(days_raw, 5));
        std::shared_ptr<arrow::Array> days;
        ARROW_ASSIGN_OR_RAISE(days, int8builder.Finish());

        int8_t months_raw[5] = {1, 3, 5, 7, 1};
        ARROW_RETURN_NOT_OK(int8builder.AppendValues(months_raw, 5));
        std::shared_ptr<arrow::Array> months;
        ARROW_ASSIGN_OR_RAISE(months, int8builder.Finish());

        arrow::Int16Builder int16builder;
        int16_t years_raw[5] = {1990, 2000, 1995, 2000, 1995};
        ARROW_RETURN_NOT_OK(int16builder.AppendValues(years_raw, 5));
        std::shared_ptr<arrow::Array> years;
        ARROW_ASSIGN_OR_RAISE(years, int16builder.Finish());

        // Get a vector of our Arrays
        std::vector<std::shared_ptr<arrow::Array>> columns = {days, months, years};

        // Make a schema to initialize the Table with
        std::shared_ptr<arrow::Field> field_day, field_month, field_year;
        std::shared_ptr<arrow::Schema> schema;

        field_day = arrow::field("Day", arrow::int8());
        field_month = arrow::field("Month", arrow::int8());
        field_year = arrow::field("Year", arrow::int16());

        schema = arrow::schema({field_day, field_month, field_year});
        // With the schema and data, create a Table
        std::shared_ptr<arrow::Table> table;
        table = arrow::Table::Make(schema, columns);

        // Write out test files in IPC, CSV, and Parquet for the example to use.
        std::shared_ptr<arrow::io::FileOutputStream> outfile;

        ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_in.arrow"));
        ARROW_ASSIGN_OR_RAISE(std::shared_ptr<arrow::ipc::RecordBatchWriter> ipc_writer,
                              arrow::ipc::MakeFileWriter(outfile, schema));

        ARROW_RETURN_NOT_OK(ipc_writer->WriteTable(*table));
        ARROW_RETURN_NOT_OK(ipc_writer->Close());

        ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_in.csv"));
        ARROW_ASSIGN_OR_RAISE(
                auto csv_writer,
                arrow::csv::MakeCSVWriter(outfile, table->schema())
        );

        ARROW_RETURN_NOT_OK(csv_writer->WriteTable(*table));
        ARROW_RETURN_NOT_OK(csv_writer->Close());

        ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_in.parquet"));
        PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 5));

        return arrow::Status::OK();
    }

}//namespace
bool MsParquetReader::checkParquetStatus() {
    const arrow::Status st = GenInitialFile();
    return st.ok();
}

namespace {

    std::map<std::string, int> buildSchemaMap(const std::shared_ptr<arrow::Schema> &schema) {

        std::string test = schema->ToString();
        const std::vector<std::string> testSplit = schema->field_names();

        std::map<std::string, int> schemaMap;

        std::cout << "PARQUET SCHEMA" << std::endl;
        for (int i = 0; i < testSplit.size(); i++) {
            std::cout << testSplit.at(i) << std::endl;
            schemaMap.emplace(testSplit.at(i), i);
        }

        return schemaMap;
    }

    template<typename T, typename U>
    bool readColumn(
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


    bool isEqual(int a, int b) {
        return a == b;
    }

    bool buildParquetRowsVector(
            std::vector<double> &mzColData,
            std::vector<float> &intensityColData,
            std::vector<int> &msLevelColData,
            std::vector<int> &scanNumberColData,
            std::vector<double> &retentionTimeColData,
            std::vector<double> &collisionEnergyColData,
            std::vector<double> &mzTargetColData,
            std::vector<double> &isoWindowLowerColData,
            std::vector<double> &isoWindowUpperColData,
            std::vector<ParquetRow> *parquetRows
            ) {

        bool e = true;

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(intensityColData.size())
                );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(msLevelColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(scanNumberColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(retentionTimeColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(collisionEnergyColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(mzTargetColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(isoWindowLowerColData.size())
        );
        if (!e) {
            return false;
        }

        e = isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(isoWindowUpperColData.size())
        );
        if (!e) {
            return false;
        }

        for (size_t i = 0; i < mzColData.size(); i++) {

            ParquetRow pr;
            pr.mz = mzColData.at(i);
            pr.intensity = intensityColData.at(i);
            pr.msLevel = msLevelColData.at(i);
            pr.scanNumber = scanNumberColData.at(i);
            pr.retentionTime = retentionTimeColData.at(i);
            pr.collisionEnergy = collisionEnergyColData.at(i);
            pr.isoWindowLower = isoWindowLowerColData.at(i);
            pr.isoWindowUpper = isoWindowUpperColData.at(i);
            pr.mzTarget = mzTargetColData.at(i);

            parquetRows->push_back(pr);
        }

        return true;
    }

}//namespace
bool MsParquetReader::readFile(const std::string &fileURI) {

    arrow::Status st;

    arrow::MemoryPool* pool = arrow::default_memory_pool();
    arrow::fs::LocalFileSystem fileSystem;

    std::shared_ptr<arrow::io::RandomAccessFile> input
            = fileSystem.OpenInputFile(fileURI).ValueOrDie();

    std::unique_ptr<parquet::arrow::FileReader> arrowReader;
    st = parquet::arrow::OpenFile(input, pool, &arrowReader);
    if (!st.ok()) {
        return false;
    }

    std::shared_ptr<arrow::Table> table;
    st = arrowReader->ReadTable(&table);
    if (!st.ok()) {
        return false;
    }

    const std::shared_ptr<arrow::Schema> schema = table->schema();
    const std::map<std::string, int> schemaMap = buildSchemaMap(schema);

    bool columnChecker = true;

    std::vector<double> mzColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("mz"), &mzColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<float> intensityColData;
    columnChecker = readColumn<arrow::FloatArray , float>(table, schema, schemaMap.at("intensity"), &intensityColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<int> msLevelColData;
    columnChecker = readColumn<arrow::Int64Array, int>(table, schema, schemaMap.at("ms_level"), &msLevelColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<int> scanNumberColData;
    columnChecker = readColumn<arrow::Int64Array, int>(table, schema, schemaMap.at("scan_nubmer"), &scanNumberColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> retentionTimeColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("retention_time"), &retentionTimeColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> collisionEnergyColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("collision_energy"), &collisionEnergyColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> mzTargetColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("mz_target"), &mzTargetColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> isoWindowLowerColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("iso_window_lower"), &isoWindowLowerColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> isoWindowUpperColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("iso_window_upper"), &isoWindowUpperColData);
    if (!columnChecker) {
        return false;
    }

    columnChecker = buildParquetRowsVector(
            mzColData,
            intensityColData,
            msLevelColData,
            scanNumberColData,
            retentionTimeColData,
            collisionEnergyColData,
            mzTargetColData,
            isoWindowLowerColData,
            isoWindowUpperColData,
            &m_parquetRows
            );

    return columnChecker;
}

bool MsParquetReader::getScans(
        int msLevel,
        QMap<ScanNumber, ScanPoints> *scanPoints
        ) {

    if (m_parquetRows.empty()) {
        return false;
    }

    scanPoints->clear();

    for (const ParquetRow &pr : m_parquetRows) {

        if (pr.msLevel != msLevel) {
            continue;
        }

        (*scanPoints)[pr.scanNumber].push_back({pr.mz, pr.intensity});
    }

    return true;
}

//Err ParquetReader::getScansInfo(std::map<ScanNumber, ScanInfo> *scansInfo) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isNotEmpty(m_parquetRows); ree;
//
//    for (const ParquetRow &pr : m_parquetRows) {
//
//        if (scansInfo->count(pr.scanNumber)) {
//            continue;
//        }
//
//        ScanInfo scanInfo;
//        scanInfo.scanNumber = pr.scanNumber;
//        scanInfo.msLevel = pr.msLevel;
//        scanInfo.retentionTime = pr.retentionTime;
//        scanInfo.collisionEnergy = pr.collisionEnergy;
//        scanInfo.mzTarget = pr.mzTarget;
//        scanInfo.isoWindowLower = pr.isoWindowLower;
//        scanInfo.isoWindowUpper = pr.isoWindowUpper;
//
//        (*scansInfo)[scanInfo.scanNumber] = scanInfo;
//    }
//
//    ERR_RETURN
//}
//
//namespace {
//
//    Err splitParquetRowsIntoSeparateVectors(
//            const std::vector<ParquetRow> &parquetRows,
//            std::vector<double> *mz,
//            std::vector<float> *intensity,
//            std::vector<int64_t> *msLevel,
//            std::vector<int64_t> *scanNumber,
//            std::vector<double> *retentionTime,
//            std::vector<double> *collisionEnergy,
//            std::vector<double> *mzTarget,
//            std::vector<double> *isoWindowLower,
//            std::vector<double> *isoWindowUpper
//    ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(parquetRows); ree;
//
//        for (const ParquetRow &pr : parquetRows) {
//            mz->push_back(pr.mz);
//            intensity->push_back(pr.intensity);
//            msLevel->push_back(static_cast<int64_t>(pr.msLevel));
//            scanNumber->push_back(static_cast<int64_t>(pr.scanNumber));
//            retentionTime->push_back(pr.retentionTime);
//            collisionEnergy->push_back(pr.collisionEnergy);
//            mzTarget->push_back(pr.mzTarget);
//            isoWindowLower->push_back(pr.isoWindowLower);
//            isoWindowUpper->push_back(pr.isoWindowUpper);
//        }
//
//        ERR_RETURN
//    }
//
//    // \brief Append a sequence of elements in one shot
//    // \Values for typenames example
//    // \arrow::Int64Builder/int64_t
//    // \arrow::DoubleBuilder/double
//    // \arrow::FloatBuilder/float
//    // arrow::StringBuilder/str
//
//    /// \brief Append a sequence of elements in one shot
//    /// \param[in] values a contiguous C array of values
//    /// \param[in] length the number of values to append
//    /// \param[in] valid_bytes an optional sequence of bytes where non-zero
//    /// indicates a valid (non-null) value
//    /// \return Status
//    template<typename ArrowBuilderType, typename T>
//    arrow::Status buildParquetDataArray(
//            const std::vector<T> &vec,
//            std::shared_ptr<arrow::Array> *output
//            ) {
//
//        arrow::Status st;
//
//        ArrowBuilderType arraybuilder;
//        ARROW_RETURN_NOT_OK(arraybuilder.AppendValues(&(vec[0]), vec.size()));
//        ARROW_ASSIGN_OR_RAISE(*output, arraybuilder.Finish());
//
//        return st;
//    }
//
//    arrow::Status writeFileLogic(
//            const std::string &outputFilePath,
//            const std::vector<ParquetRow> &parquetRows
//            ) {
//
//        std::vector<double> mz;
//        std::vector<float> intensity;
//        std::vector<int64_t> msLevel;
//        std::vector<int64_t> scanNumber;
//        std::vector<double> retentionTime;
//        std::vector<double> collisionEnergy;
//        std::vector<double> mzTarget;
//        std::vector<double> isoWindowLower;
//        std::vector<double> isoWindowUpper;
//
//        arrow::Status st;
//
//        Err e = splitParquetRowsIntoSeparateVectors(
//                parquetRows,
//                &mz,
//                &intensity,
//                &msLevel,
//                &scanNumber,
//                &retentionTime,
//                &collisionEnergy,
//                &mzTarget,
//                &isoWindowLower,
//                &isoWindowUpper
//                );
//
//        if (e != Error::eNoError) {
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> mzArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(mz, &mzArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//
//        std::shared_ptr<arrow::Array> intensityArrow;
//        st = buildParquetDataArray<arrow::FloatBuilder>(intensity, &intensityArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> msLevelArrow;
//        st = buildParquetDataArray<arrow::Int64Builder>(msLevel, &msLevelArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> scanNumberArrow;
//        st = buildParquetDataArray<arrow::Int64Builder>(scanNumber, &scanNumberArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> retentionTimeArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(retentionTime, &retentionTimeArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> collisionEnergyArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(collisionEnergy, &collisionEnergyArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> mzTargetArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(mzTarget, &mzTargetArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> isoWindowLowerArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(isoWindowLower, &isoWindowLowerArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> isoWindowUpperArrow;
//        st = buildParquetDataArray<arrow::DoubleBuilder>(isoWindowUpper, &isoWindowUpperArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::vector<std::shared_ptr<arrow::Array>> columns = {
//                mzArrow,
//                intensityArrow,
//                msLevelArrow,
//                scanNumberArrow,
//                retentionTimeArrow,
//                collisionEnergyArrow,
//                mzTargetArrow,
//                isoWindowLowerArrow,
//                isoWindowUpperArrow
//        };
//
//        std::shared_ptr<arrow::Field> mzField;
//        std::shared_ptr<arrow::Field> intensityField;
//        std::shared_ptr<arrow::Field> msLevelField;
//        std::shared_ptr<arrow::Field> scanNumberField;
//        std::shared_ptr<arrow::Field> retentionTimeField;
//        std::shared_ptr<arrow::Field> collisionEnergyField;
//        std::shared_ptr<arrow::Field> mzTargetField;
//        std::shared_ptr<arrow::Field> isoWindowLowerField;
//        std::shared_ptr<arrow::Field> isoWindowUpperField;
//
//        std::shared_ptr<arrow::Schema> schema;
//
//        mzField = arrow::field("mz", arrow::float64());
//        intensityField = arrow::field("intensity", arrow::float32());
//        msLevelField = arrow::field("ms_level", arrow::int64());
//        scanNumberField = arrow::field("scan_number", arrow::int64());
//        retentionTimeField = arrow::field("retention_time", arrow::float64());
//        collisionEnergyField = arrow::field("collision_energy", arrow::float64());
//        mzTargetField = arrow::field("mz_target", arrow::float64());
//        isoWindowLowerField = arrow::field("iso_window_lower", arrow::float64());
//        isoWindowUpperField = arrow::field("iso_window_upper", arrow::float64());
//
//        schema = arrow::schema({
//            mzField,
//            intensityField,
//            msLevelField,
//            scanNumberField,
//            retentionTimeField,
//            collisionEnergyField,
//            mzTargetField,
//            isoWindowLowerField,
//            isoWindowUpperField
//        });
//
//        std::shared_ptr<arrow::Table> table;
//        table = arrow::Table::Make(schema, columns);
//
//        std::shared_ptr<arrow::io::FileOutputStream> outfile;
//
//        ARROW_ASSIGN_OR_RAISE(
//                outfile,
//                arrow::io::FileOutputStream::Open(outputFilePath)
//                );
//
//        PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
//                *table,
//                arrow::default_memory_pool(),
//                outfile
//                ));
//
//        return st;
//    }
//
//}//namespace
//Err ParquetReader::writeFile(const std::string &outputFilePath){
//
//    ERR_INIT
//
//    arrow::Status st =  writeFileLogic(outputFilePath, m_parquetRows);
//    e = ErrorUtils::isTrue(st.ok()); ree;
//
//    ERR_RETURN
//}
//
//void ParquetReader::clearParquetRows() {
//    m_parquetRows.clear();
//}
//
//void ParquetReader::addParquetRowToParquetRows(const ParquetRow &parquetRow) {
//    m_parquetRows.push_back(parquetRow);
//}
