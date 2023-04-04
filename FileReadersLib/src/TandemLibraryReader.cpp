//
// Created by anichols on 4/3/23.
//

#include "TandemLibraryReader.h"

#include "ErrorUtils.h"


namespace {

    bool splitParquetRowsIntoSeparateVectors(
            const std::vector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
            std::vector<std::string> *peptideSequences,
            std::vector<std::string> *intensityVecsByteString,
            std::vector<std::string> *ionLabelVecs
    ) {

//        if (tandemLibraryReaderRows.empty()) {
//            return false;
//        }
//
//        for (const TandemLibraryReaderRow &tlrs : tandemLibraryReaderRows) {
//            peptideSequences->push_back(tlrs.peptideString.toStdString());
//            intensityVecsByteString->push_back(tlrs.intensityVals);
//            ionLabelVecs->push_back(tlrs.ionLabels);
//        }

        return true;
    }

    arrow::Status writeFileLogic(
            const std::string &outputFilePath,
            const std::vector<TandemLibraryReaderRow> &tandemLibraryReaderRows
    ) {

        std::vector<std::string> peptideSequences;
        std::vector<std::vector<double>> intensityVecs;
        std::vector<std::vector<std::string>> ionLabelVecs;

        arrow::Status st;

//        bool e = splitParquetRowsIntoSeparateVectors(
//                tandemLibraryReaderRows,
//                &peptideSequences,
//                &intensityVecs,
//                &ionLabelVecs
//        );

//        if (!e) {
//            st;
//        }
//
//        std::shared_ptr<arrow::Array> mzArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::Str>(mz, &mzArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> intensityArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::FloatBuilder>(intensity, &intensityArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> msLevelArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::Int64Builder>(msLevel, &msLevelArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> scanNumberArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::Int64Builder>(scanNumber, &scanNumberArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> retentionTimeArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(retentionTime, &retentionTimeArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> collisionEnergyArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(collisionEnergy, &collisionEnergyArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> mzTargetArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(mzTarget, &mzTargetArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> isoWindowLowerArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(isoWindowLower, &isoWindowLowerArrow);
//        if (!st.ok()){
//            return st;
//        }
//
//        std::shared_ptr<arrow::Array> isoWindowUpperArrow;
//        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(isoWindowUpper, &isoWindowUpperArrow);
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
//                                       mzField,
//                                       intensityField,
//                                       msLevelField,
//                                       scanNumberField,
//                                       retentionTimeField,
//                                       collisionEnergyField,
//                                       mzTargetField,
//                                       isoWindowLowerField,
//                                       isoWindowUpperField
//                               });
//
//        std::shared_ptr<arrow::Table> table;
//        table = arrow::Table::Make(schema, columns);
//
//        std::shared_ptr<arrow::io::FileOutputStream> outfile;
//
//        ARROW_ASSIGN_OR_RAISE(
//                outfile,
//                arrow::io::FileOutputStream::Open(outputFilePath)
//        );
//
//        PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
//                *table,
//                arrow::default_memory_pool(),
//                outfile
//        ));

        return st;
    }



}//namespace
Err TandemLibraryReader::writeTandemPrediction(
        const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
        const QString &outputFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(tandemLibraryReaderRows); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;




    ERR_RETURN
}
