//
// Created by anichols on 4/3/23.
//

#include "TandemLibraryReader.h"

#include "ErrorUtils.h"


namespace {

    const std::string PEPTIDE_SEQUENCE = "peptideSequence";
    const std::string INTENSITY_VEC = "intensityVec";
    const std::string ION_LABELS = "ionLabels";

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
            return st;
        }

        std::shared_ptr<arrow::Array> peptideSequencesArrow;
        st = ParquetReaderBase::buildParquetDataArrayFromString(peptideSequences, &peptideSequencesArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> intensityVecsArrow;
        st = ParquetReaderBase::buildParquetDataArrayFromByteArray(intensityVecs, &intensityVecsArrow);
        if (!st.ok()){
            return st;
        }

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

        peptideSequencesField = arrow::field(PEPTIDE_SEQUENCE, arrow::utf8());
        intensityVecsField = arrow::field(INTENSITY_VEC, arrow::binary());
        ionLabelVecsField = arrow::field(ION_LABELS, arrow::utf8());

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

namespace {

    Err buildTandemLibraryReaderRowVector(
            std::vector<std::string> &peptideSequenceColData,
            std::vector<std::string> &intensityColData,
            std::vector<std::string> &ionLabelsColData,
            QVector<TandemLibraryReaderRow> *tandemLibraryReaderRows
    ) {

        ERR_INIT

        bool ee = true;
        ee = ParquetReaderBase::isEqual(
                static_cast<int>(peptideSequenceColData.size()),
                static_cast<int>(intensityColData.size())
        );
        if (!ee) {
            rrr(eError);
        }

        ee = ParquetReaderBase::isEqual(
                static_cast<int>(peptideSequenceColData.size()),
                static_cast<int>(ionLabelsColData.size())
        );
        if (!ee) {
            return eError;
        }

        for (size_t i = 0; i < peptideSequenceColData.size(); i++) {

            TandemLibraryReaderRow tlrr;
            tlrr.peptideString = QString::fromStdString(peptideSequenceColData.at(i));
            tlrr.intensityVals =  ParquetReaderBase::bytesStdStringToQVector<double>(intensityColData.at(i));
            tlrr.ionLabels = QString::fromStdString(ionLabelsColData.at(i)).split(S_GLOBAL_SETTINGS.SEPARATOR);

            tandemLibraryReaderRows->push_back(tlrr);
        }

        ERR_RETURN
    }

}//namespace
Err TandemLibraryReader::readTandemPredictions(
        const std::string &fileURI,
        QVector<TandemLibraryReaderRow> *tandemLibraryReaderRows
        ) {

    ERR_INIT

    arrow::Status st;

    arrow::MemoryPool* pool = arrow::default_memory_pool();
    arrow::fs::LocalFileSystem fileSystem;

    std::shared_ptr<arrow::io::RandomAccessFile> input
            = fileSystem.OpenInputFile(fileURI).ValueOrDie();

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

    const std::shared_ptr<arrow::Schema> schema = table->schema();
    const std::map<std::string, int> schemaMap = ParquetReaderBase::buildSchemaMap(schema);

    bool columnChecker = true;

    std::vector<std::string> peptideSequenceColData;
    columnChecker = readColumn<arrow::StringArray , std::string>(
            table,
            schema,
            schemaMap.at(PEPTIDE_SEQUENCE),
            &peptideSequenceColData
            );
    if (!columnChecker) {
        rrr(eError);
    }

    std::vector<std::string> intensityColData;
    columnChecker = readColumn<arrow::BinaryArray, std::string>(
            table,
            schema,
            schemaMap.at(INTENSITY_VEC),
            &intensityColData
            );
    if (!columnChecker) {
        rrr(eError);
    }

    std::vector<std::string> ionLabelsColData;
    columnChecker = readColumn<arrow::StringArray, std::string>(
            table,
            schema,
            schemaMap.at(ION_LABELS),
            &ionLabelsColData
            );
    if (!columnChecker) {
        rrr(eError);
    }

    e = buildTandemLibraryReaderRowVector(
            peptideSequenceColData,
            intensityColData,
            ionLabelsColData,
            tandemLibraryReaderRows
    );

    ERR_RETURN
}
