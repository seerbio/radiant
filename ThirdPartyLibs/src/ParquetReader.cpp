//
// Created by anichols on 2/6/23.
//

#include "ParquetReader.h"

#include "ErrorUtils.h"
#include "StringUtils.h"

#include <arrow/api.h>
#include <arrow/csv/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>


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
Err ParquetReader::checkParquetStatus() {

    ERR_INIT

    const arrow::Status st = GenInitialFile();
    e = ErrorUtils::isTrue(st.ok()); ree;

    ERR_RETURN
}

namespace {

    void printSchema(const std::shared_ptr<arrow::Schema> &schema) {

        std::string test = schema->ToString();
        const std::vector<std::string> testSplit = schema->field_names();

        std::cout << "PARQUET SCHEMA" << std::endl;
        for (auto &s : testSplit) {
            std::cout << s << std::endl;
        }
    }

    template<typename T, typename U>
    Err readColumn(
            const std::shared_ptr<arrow::Table> &table,
            const std::shared_ptr<arrow::Schema> &schema,
            int columnIndex,
            std::vector<U> *output
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(schema->field_names()); ree;

        std::shared_ptr<arrow::ChunkedArray> mzCol = table->column(columnIndex);
        const std::shared_ptr<arrow::Array> mzColChunks = mzCol->chunks()[0];

        auto arrowArray = std::static_pointer_cast<T>(mzColChunks);

        for (int64_t i = 0; i < mzColChunks->length(); ++i) {
            output->push_back(arrowArray->Value(i));
        }

        ERR_RETURN
    }


    Err buildParquetRowsVector(
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

        ERR_INIT

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(intensityColData.size())
                ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(msLevelColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(scanNumberColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(retentionTimeColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(collisionEnergyColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(mzTargetColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(isoWindowLowerColData.size())
        ); ree;

        e = ErrorUtils::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(isoWindowUpperColData.size())
        ); ree;

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

        ERR_RETURN
    }

}//namespace
Err ParquetReader::readFile(const std::string &fileURI) {

    ERR_INIT

    e = checkParquetStatus(); ree;

    arrow::Status st;

    arrow::MemoryPool* pool = arrow::default_memory_pool();
    arrow::fs::LocalFileSystem fileSystem;

    std::shared_ptr<arrow::io::RandomAccessFile> input
            = fileSystem.OpenInputFile(fileURI).ValueOrDie();

    std::unique_ptr<parquet::arrow::FileReader> arrowReader;
    st = parquet::arrow::OpenFile(input, pool, &arrowReader);
    e = ErrorUtils::isTrue(st.ok()); ree;

    std::shared_ptr<arrow::Table> table;
    st = arrowReader->ReadTable(&table);
    e = ErrorUtils::isTrue(st.ok()); ree;

    const std::shared_ptr<arrow::Schema> schema = table->schema();
//    printSchema(schema);

    std::vector<double> mzColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 0, &mzColData);

    std::vector<float> intensityColData;
    e = readColumn<arrow::FloatArray , float>(table, schema, 1, &intensityColData);

    std::vector<int> msLevelColData;
    e = readColumn<arrow::Int64Array, int>(table, schema, 2, &msLevelColData);

    std::vector<int> scanNumberColData;
    e = readColumn<arrow::Int64Array, int>(table, schema, 3, &scanNumberColData);

    std::vector<double> retentionTimeColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 4, &retentionTimeColData);

    std::vector<double> collisionEnergyColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 5, &collisionEnergyColData);

    std::vector<double> mzTargetColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 6, &mzTargetColData);

    std::vector<double> isoWindowLowerColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 7, &isoWindowLowerColData);

    std::vector<double> isoWindowUpperColData;
    e = readColumn<arrow::DoubleArray, double>(table, schema, 8, &isoWindowUpperColData);

    e = buildParquetRowsVector(
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
            ); ree;

    ERR_RETURN
}

Err ParquetReader::getScans(
        int msLevel,
        std::map<ScanNumber, ScanPoints> *scanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_parquetRows); ree;

    scanPoints->clear();

    for (const ParquetRow &pr : m_parquetRows) {

        if (pr.msLevel != msLevel) {
            continue;
        }

        (*scanPoints)[pr.scanNumber].emplace_back(pr.mz, pr.intensity);
    }

    ERR_RETURN
}

Err ParquetReader::getScansInfo(std::map<ScanNumber, ScanInfo> *scansInfo) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_parquetRows); ree;

    for (const ParquetRow &pr : m_parquetRows) {

        if (scansInfo->count(pr.scanNumber)) {
            continue;
        }

        ScanInfo scanInfo;
        scanInfo.scanNumber = pr.scanNumber;
        scanInfo.msLevel = pr.msLevel;
        scanInfo.retentionTime = pr.retentionTime;
        scanInfo.collisionEnergy = pr.collisionEnergy;
        scanInfo.mzTarget = pr.mzTarget;
        scanInfo.isoWindowLower = pr.isoWindowLower;
        scanInfo.isoWindowUpper = pr.isoWindowUpper;

        (*scansInfo)[scanInfo.scanNumber] = scanInfo;
    }

    ERR_RETURN
}

namespace {

    Err splitParquetRowsIntoSeparateVectors(
            const std::vector<ParquetRow> &parquetRows,
            std::vector<double> *mz,
            std::vector<float> *intensity,
            std::vector<int64_t> *msLevel,
            std::vector<int64_t> *scanNumber,
            std::vector<double> *retentionTime,
            std::vector<double> *collisionEnergy,
            std::vector<double> *mzTarget,
            std::vector<double> *isoWindowLower,
            std::vector<double> *isoWindowUpper
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(parquetRows); ree;

        for (const ParquetRow &pr : parquetRows) {
            mz->push_back(pr.mz);
            intensity->push_back(pr.intensity);
            msLevel->push_back(static_cast<int64_t>(pr.msLevel));
            scanNumber->push_back(static_cast<int64_t>(pr.scanNumber));
            retentionTime->push_back(pr.retentionTime);
            collisionEnergy->push_back(pr.collisionEnergy);
            mzTarget->push_back(pr.mzTarget);
            isoWindowLower->push_back(pr.isoWindowLower);
            isoWindowUpper->push_back(pr.isoWindowUpper);
        }

        ERR_RETURN
    }

    // \brief Append a sequence of elements in one shot
    // \Values for typenames example
    // \arrow::Int64Builder/int64_t
    // \arrow::DoubleBuilder/double
    // \arrow::FloatBuilder/float
    // arrow::StringBuilder/str

    /// \brief Append a sequence of elements in one shot
    /// \param[in] values a contiguous C array of values
    /// \param[in] length the number of values to append
    /// \param[in] valid_bytes an optional sequence of bytes where non-zero
    /// indicates a valid (non-null) value
    /// \return Status
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

    arrow::Status writeFileLogic(
            const std::string &outputFilePath,
            const std::vector<ParquetRow> &parquetRows
            ) {

        std::vector<double> mz;
        std::vector<float> intensity;
        std::vector<int64_t> msLevel;
        std::vector<int64_t> scanNumber;
        std::vector<double> retentionTime;
        std::vector<double> collisionEnergy;
        std::vector<double> mzTarget;
        std::vector<double> isoWindowLower;
        std::vector<double> isoWindowUpper;

        arrow::Status st;

        Err e = splitParquetRowsIntoSeparateVectors(
                parquetRows,
                &mz,
                &intensity,
                &msLevel,
                &scanNumber,
                &retentionTime,
                &collisionEnergy,
                &mzTarget,
                &isoWindowLower,
                &isoWindowUpper
                );

        if (e != Error::eNoError) {
            return st;
        }

        std::shared_ptr<arrow::Array> mzArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(mz, &mzArrow);
        if (!st.ok()){
            return st;
        }


        std::shared_ptr<arrow::Array> intensityArrow;
        st = buildParquetDataArray<arrow::FloatBuilder>(intensity, &intensityArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> msLevelArrow;
        st = buildParquetDataArray<arrow::Int64Builder>(msLevel, &msLevelArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> scanNumberArrow;
        st = buildParquetDataArray<arrow::Int64Builder>(scanNumber, &scanNumberArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> retentionTimeArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(retentionTime, &retentionTimeArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> collisionEnergyArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(collisionEnergy, &collisionEnergyArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> mzTargetArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(mzTarget, &mzTargetArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> isoWindowLowerArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(isoWindowLower, &isoWindowLowerArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> isoWindowUpperArrow;
        st = buildParquetDataArray<arrow::DoubleBuilder>(isoWindowUpper, &isoWindowUpperArrow);
        if (!st.ok()){
            return st;
        }

        std::vector<std::shared_ptr<arrow::Array>> columns = {
                mzArrow,
                intensityArrow,
                msLevelArrow,
                scanNumberArrow,
                retentionTimeArrow,
                collisionEnergyArrow,
                mzTargetArrow,
                isoWindowLowerArrow,
                isoWindowUpperArrow
        };

        std::shared_ptr<arrow::Field> mzField;
        std::shared_ptr<arrow::Field> intensityField;
        std::shared_ptr<arrow::Field> msLevelField;
        std::shared_ptr<arrow::Field> scanNumberField;
        std::shared_ptr<arrow::Field> retentionTimeField;
        std::shared_ptr<arrow::Field> collisionEnergyField;
        std::shared_ptr<arrow::Field> mzTargetField;
        std::shared_ptr<arrow::Field> isoWindowLowerField;
        std::shared_ptr<arrow::Field> isoWindowUpperField;

        std::shared_ptr<arrow::Schema> schema;

        mzField = arrow::field("mz", arrow::float64());
        intensityField = arrow::field("intensity", arrow::float32());
        msLevelField = arrow::field("ms_level", arrow::int64());
        scanNumberField = arrow::field("scan_number", arrow::int64());
        retentionTimeField = arrow::field("retention_time", arrow::float64());
        collisionEnergyField = arrow::field("collision_energy", arrow::float64());
        mzTargetField = arrow::field("mz_target", arrow::float64());
        isoWindowLowerField = arrow::field("iso_window_lower", arrow::float64());
        isoWindowUpperField = arrow::field("iso_window_upper", arrow::float64());

        schema = arrow::schema({
            mzField,
            intensityField,
            msLevelField,
            scanNumberField,
            retentionTimeField,
            collisionEnergyField,
            mzTargetField,
            isoWindowLowerField,
            isoWindowUpperField
        });

        std::shared_ptr<arrow::Table> table;
        table = arrow::Table::Make(schema, columns);

        std::shared_ptr<arrow::io::FileOutputStream> outfile;

        ARROW_ASSIGN_OR_RAISE(
                outfile,
                arrow::io::FileOutputStream::Open(outputFilePath)
                );

        PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
                *table,
                arrow::default_memory_pool(),
                outfile
                ));

        return st;
    }

}//namespace
Err ParquetReader::writeFile(const std::string &outputFilePath){

    ERR_INIT

    arrow::Status st =  writeFileLogic(outputFilePath, m_parquetRows);
    e = ErrorUtils::isTrue(st.ok()); ree;

    ERR_RETURN
}

void ParquetReader::clearParquetRows() {
    m_parquetRows.clear();
}

void ParquetReader::addParquetRowToParquetRows(const ParquetRow &parquetRow) {
    m_parquetRows.push_back(parquetRow);
}
