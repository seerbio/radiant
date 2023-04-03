//
// Created by anichols on 2/6/23.
//

#include "MsParquetReader.h"

//#include "ErrorUtils.h"
//#include "StringUtils.h"
#include <QVector>
#include <QPointF>

#include <iostream>


namespace {

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

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(intensityColData.size())
                );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(msLevelColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(scanNumberColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(retentionTimeColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(collisionEnergyColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(mzTargetColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(mzColData.size()),
                static_cast<int>(isoWindowLowerColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
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
    const std::map<std::string, int> schemaMap = ParquetReaderBase::buildSchemaMap(schema);

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

bool MsParquetReader::getScansInfo(QMap<ScanNumber, ScanInfo> *scansInfo) {

    if (m_parquetRows.empty()) {
        return false;
    }

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

    return true;
}

namespace {

    bool splitParquetRowsIntoSeparateVectors(
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

        if (parquetRows.empty()) {
            return false;
        }

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

        return true;
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

        bool e = splitParquetRowsIntoSeparateVectors(
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

        if (!e) {
            st;
        }

        std::shared_ptr<arrow::Array> mzArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(mz, &mzArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> intensityArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::FloatBuilder>(intensity, &intensityArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> msLevelArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::Int64Builder>(msLevel, &msLevelArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> scanNumberArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::Int64Builder>(scanNumber, &scanNumberArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> retentionTimeArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(retentionTime, &retentionTimeArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> collisionEnergyArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(collisionEnergy, &collisionEnergyArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> mzTargetArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(mzTarget, &mzTargetArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> isoWindowLowerArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(isoWindowLower, &isoWindowLowerArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> isoWindowUpperArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(isoWindowUpper, &isoWindowUpperArrow);
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
bool MsParquetReader::writeFile(const std::string &outputFilePath){

    arrow::Status st =  writeFileLogic(outputFilePath, m_parquetRows);
    if (!st.ok()) {
        return false;
    }

    return true;
}

void MsParquetReader::clearParquetRows() {
    m_parquetRows.clear();
}

void MsParquetReader::addParquetRowToParquetRows(const ParquetRow &parquetRow) {
    m_parquetRows.push_back(parquetRow);
}
