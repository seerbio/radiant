//
// Created by anichols on 4/3/23.
//

#include "TandemLibraryReader.h"

#include "ErrorUtils.h"


namespace {

    bool splitParquetRowsIntoSeparateVectors(
            const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
            std::vector<std::string> *peptideSequences,
            std::vector<std::string> *intensityVecsByteString,
            std::vector<std::string> *ionLabelVecs
    ) {

        if (tandemLibraryReaderRows.empty()) {
            return false;
        }

        for (const TandemLibraryReaderRow &tlrs : tandemLibraryReaderRows) {

            const std::string intensityByteStream
                    = ParquetReaderBase::qVectorToBytesStdString(tlrs.intensityVals);

            const std::string ionLabelsJoined
                    = tlrs.ionLabels.join(S_GLOBAL_SETTINGS.SEPARATOR).toStdString();

            peptideSequences->push_back(tlrs.peptideString.toStdString());
            intensityVecsByteString->push_back(intensityByteStream);
            ionLabelVecs->push_back(ionLabelsJoined);
        }

        return true;
    }

    arrow::Status writeFileLogic(
            const std::string &outputFilePath,
            const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows
    ) {

        std::vector<std::string> peptideSequences;
        std::vector<std::string> intensityVecs;
        std::vector<std::string> ionLabelVecs;

        arrow::Status st;

        bool e = splitParquetRowsIntoSeparateVectors(
                tandemLibraryReaderRows,
                &peptideSequences,
                &intensityVecs,
                &ionLabelVecs
        );


        if (!e) {
            qDebug() << "SDFKSDF" << st.ok();
            return st;
        }

        std::shared_ptr<arrow::Array> peptideSequencesArrow;
        st = ParquetReaderBase::buildParquetDataArrayFromString(peptideSequences, &peptideSequencesArrow);
        if (!st.ok()){
            qDebug() << "KDKFDJ" << st.ok();
            return st;
        }

        std::shared_ptr<arrow::Array> intensityVecsArrow;
        st = ParquetReaderBase::buildParquetDataArrayFromByteArray(intensityVecs, &intensityVecsArrow);
        if (!st.ok()){
            return st;
        }

//        std::shared_ptr<arrow::Array> intensityVecsArrow;
//        st = ParquetReaderBase::buildParquetDataArrayFromString(intensityVecs, &intensityVecsArrow);
//        if (!st.ok()){
//            return st;
//        }

        std::shared_ptr<arrow::Array> ionLabelVecsArrow;
        st = ParquetReaderBase::buildParquetDataArrayFromString(ionLabelVecs, &ionLabelVecsArrow);
        if (!st.ok()){
            return st;
        }


        std::vector<std::shared_ptr<arrow::Array>> columns = {
                peptideSequencesArrow,
                intensityVecsArrow,
                ionLabelVecsArrow
        };

        std::shared_ptr<arrow::Field> peptideSequencesField;
        std::shared_ptr<arrow::Field> intensityVecsField;
        std::shared_ptr<arrow::Field> ionLabelVecsField;

        std::shared_ptr<arrow::Schema> schema;

        peptideSequencesField = arrow::field("peptideSequence", arrow::utf8());
        intensityVecsField = arrow::field("intensityVec", arrow::binary());
        ionLabelVecsField = arrow::field("ionLabels", arrow::utf8());

        schema = arrow::schema({
                    peptideSequencesField,
                    intensityVecsField,
                    ionLabelVecsField
                });

        std::shared_ptr<arrow::Table> table;
        table = arrow::Table::Make(schema, columns);

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
Err TandemLibraryReader::writeTandemPredictions(
        const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
        const QString &outputFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(tandemLibraryReaderRows); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;

    arrow::Status st = writeFileLogic(
            outputFilePath.toStdString(),
            tandemLibraryReaderRows
            );

    e = ErrorUtils::isTrue(st.ok()); ree;
    qDebug() << "File written to:" << outputFilePath;

    ERR_RETURN
}
