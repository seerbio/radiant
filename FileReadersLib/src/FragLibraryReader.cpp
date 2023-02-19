//
// Created by anichols on 2/18/23.
//

#include "FragLibraryReader.h"


namespace {

    bool buildFragIonsVector(
            std::vector<int> &peptideIdColData,
            std::vector<double> &mzFragColData,
            std::vector<double> &peptideMassColData,
            std::vector<FragLibIon> *fragLibIons
    ) {

        bool e = true;

        e = ParquetReaderBase::isEqual(
                static_cast<int>(peptideIdColData.size()),
                static_cast<int>(mzFragColData.size())
        );
        if (!e) {
            return false;
        }

        e = ParquetReaderBase::isEqual(
                static_cast<int>(peptideIdColData.size()),
                static_cast<int>(peptideMassColData.size())
        );
        if (!e) {
            return false;
        }

        for (size_t i = 0; i < peptideIdColData.size(); i++) {

            FragLibIon fli;
            fli.peptideId = peptideIdColData.at(i);
            fli.mzFrag = mzFragColData.at(i);
            fli.peptideMass = peptideMassColData.at(i);

            fragLibIons->push_back(fli);
        }

        return true;
    }

}//namespace
bool FragLibraryReader::readFile(const std::string &fileURI) {
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

    std::vector<int> peptideIdColData;
    columnChecker = readColumn<arrow::Int64Array, int>(table, schema, schemaMap.at("peptideId"), &peptideIdColData);
    if (!columnChecker) {
        return false;
    }
    qDebug() << peptideIdColData;

    std::vector<double> mzFragColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("mzFrag"), &mzFragColData);
    if (!columnChecker) {
        return false;
    }

    std::vector<double> peptideMassColData;
    columnChecker = readColumn<arrow::DoubleArray, double>(table, schema, schemaMap.at("peptideMass"), &peptideMassColData);
    if (!columnChecker) {
        return false;
    }

    columnChecker = buildFragIonsVector(
            peptideIdColData,
            mzFragColData,
            peptideMassColData,
            &m_fragLibIons
    );

    return columnChecker;
}


namespace {

    bool splitParquetRowsIntoSeparateVectors(
            const std::vector<FragLibIon> &fragLibIons,
            std::vector<int64_t> *peptideId,
            std::vector<double> *mzFrag,
            std::vector<double> *peptideMass
            ) {

        if (fragLibIons.empty()) {
            return false;
        }

        for (const FragLibIon &fli : fragLibIons) {
            peptideId->push_back(static_cast<int64_t>(fli.peptideId));
            mzFrag->push_back(fli.mzFrag);
            peptideMass->push_back(fli.peptideMass);
        }

        return true;
    }

    arrow::Status writeFileLogic(
            const std::string &outputFilePath,
            const std::vector<FragLibIon> &fragLibIons
            ) {

        std::vector<int64_t> peptideId;
        std::vector<double> mzFrag;
        std::vector<double> peptideMass;

        arrow::Status st;

        bool e = splitParquetRowsIntoSeparateVectors(
                fragLibIons,
                &peptideId,
                &mzFrag,
                &peptideMass
                );

        if (!e) {
            st;
        }

        std::shared_ptr<arrow::Array> peptideIdArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::Int64Builder>(peptideId, &peptideIdArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> mzFragArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(mzFrag, &mzFragArrow);
        if (!st.ok()){
            return st;
        }

        std::shared_ptr<arrow::Array> peptideMassArrow;
        st = ParquetReaderBase::buildParquetDataArray<arrow::DoubleBuilder>(peptideMass, &peptideMassArrow);
        if (!st.ok()){
            return st;
        }

        std::vector<std::shared_ptr<arrow::Array>> columns = {
                peptideIdArrow,
                mzFragArrow,
                peptideMassArrow
        };

        std::shared_ptr<arrow::Field> peptideIdField;
        std::shared_ptr<arrow::Field> mzFragField;
        std::shared_ptr<arrow::Field> peptideMassField;

        std::shared_ptr<arrow::Schema> schema;

        peptideIdField = arrow::field("peptideId", arrow::int64());
        mzFragField = arrow::field("mzFrag", arrow::float64());
        peptideMassField = arrow::field("peptideMass", arrow::float64());

        schema = arrow::schema({
                               peptideIdField,
                               mzFragField,
                               peptideMassField
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
bool FragLibraryReader::writeFile(const std::string &outputFilePath) {
    arrow::Status st =  writeFileLogic(outputFilePath, m_fragLibIons);
    if (!st.ok()) {
        return false;
    }

    return true;
}
