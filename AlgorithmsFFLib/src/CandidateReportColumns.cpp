// Arrow headers must precede Qt's macros.
#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include "CandidateReportColumns.h"

#include <QDir>
#include <QFile>
#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>

namespace {
template<typename T>
std::shared_ptr<arrow::DataType> fieldType() {
    if constexpr (std::is_same_v<T, float>) return arrow::float32();
    else if constexpr (std::is_same_v<T, double>) return arrow::float64();
    else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, bool>) return arrow::int64();
    else {
        static_assert(std::is_base_of_v<QString, T>, "Unsupported candidate report field type");
        return arrow::utf8();
    }
}

struct Field {
    QString name;
    std::shared_ptr<arrow::DataType> type;
    bool boolean = false;
};

std::vector<Field> reportFields() {
    std::vector<Field> fields;
    CandidateScoresReaderRow().visitFields([&](const QString &name, const auto &value) {
        using T = std::decay_t<decltype(value)>;
        fields.push_back({name, fieldType<T>(), std::is_same_v<T, bool>});
    });
    return fields;
}

arrow::Status writeTable(const arrow::Table &table, const QString &path) {
    ARROW_ASSIGN_OR_RAISE(auto output, arrow::io::FileOutputStream::Open(path.toStdString()));
    ARROW_RETURN_NOT_OK(parquet::arrow::WriteTable(table, arrow::default_memory_pool(), output));
    return output->Close();
}

arrow::Status writeRows(const QVector<CandidateScoresReaderRow> &rows, const QString &path) {
    auto fields = reportFields();
    std::vector<std::unique_ptr<arrow::ArrayBuilder>> builders;
    for (const auto &field : fields) {
        std::unique_ptr<arrow::ArrayBuilder> builder;
        ARROW_RETURN_NOT_OK(arrow::MakeBuilder(arrow::default_memory_pool(), field.type, &builder));
        ARROW_RETURN_NOT_OK(builder->Reserve(rows.size()));
        builders.push_back(std::move(builder));
    }
    for (const auto &row : rows) {
        size_t column = 0;
        arrow::Status status = arrow::Status::OK();
        row.visitFields([&](const QString &, const auto &value) {
            auto *builder = builders[column++].get();
            if (!status.ok()) return;
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, float>)
                static_cast<arrow::FloatBuilder*>(builder)->UnsafeAppend(value);
            else if constexpr (std::is_same_v<T, double>)
                static_cast<arrow::DoubleBuilder*>(builder)->UnsafeAppend(value);
            else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, bool>)
                static_cast<arrow::Int64Builder*>(builder)->UnsafeAppend(int64_t(value));
            else {
                const auto bytes = value.toUtf8();
                status = static_cast<arrow::StringBuilder*>(builder)->Append(
                    bytes.constData(), bytes.size());
            }
        });
        ARROW_RETURN_NOT_OK(status);
    }
    std::vector<size_t> order(fields.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return fields[a].name < fields[b].name; });
    std::vector<std::shared_ptr<arrow::Field>> schema;
    std::vector<std::shared_ptr<arrow::Array>> columns;
    for (size_t index : order) {
        schema.push_back(arrow::field(fields[index].name.toStdString(), fields[index].type));
        std::shared_ptr<arrow::Array> column;
        ARROW_RETURN_NOT_OK(builders[index]->Finish(&column));
        columns.push_back(std::move(column));
    }
    return writeTable(*arrow::Table::Make(arrow::schema(schema), columns), path);
}

arrow::Result<std::shared_ptr<arrow::Table>> readTable(const QString &path) {
    ARROW_ASSIGN_OR_RAISE(auto input, arrow::io::ReadableFile::Open(path.toStdString()));
    ARROW_ASSIGN_OR_RAISE(auto reader, parquet::arrow::OpenFile(input, arrow::default_memory_pool()));
    std::shared_ptr<arrow::Table> table;
    ARROW_RETURN_NOT_OK(reader->ReadTable(&table));
    ARROW_RETURN_NOT_OK(table->ValidateFull());
    return table->CombineChunks();
}

bool canonicalTable(const arrow::Table &table, const std::vector<Field> &fields) {
    if (table.num_columns() != int(fields.size())) return false;
    for (int index = 0; index < table.num_columns(); ++index) {
        const auto &field = fields[index];
        const auto &column = table.column(index);
        const auto expected = arrow::field(field.name.toStdString(), field.type);
        if (!table.field(index)->Equals(expected, false)
            || column->null_count() != 0 || column->num_chunks() != 1) return false;
        // The legacy row conversion narrows integers and canonicalizes bools.
        // Use it for noncanonical values even if the Arrow types match.
        if (field.type->id() == arrow::Type::INT64) {
            const auto values = std::static_pointer_cast<arrow::Int64Array>(column->chunk(0));
            for (int64_t row = 0; row < values->length(); ++row) {
                const auto value = values->Value(row);
                if (value < std::numeric_limits<int>::min()
                    || value > std::numeric_limits<int>::max()
                    || (field.boolean && value != 0 && value != 1)) return false;
            }
        }
    }
    return true;
}

template<typename Array>
std::shared_ptr<Array> column(const arrow::Table &table, const QString &name) {
    return std::static_pointer_cast<Array>(table.GetColumnByName(name.toStdString())->chunk(0));
}

QString text(const arrow::StringArray &array, int64_t row) {
    const auto bytes = array.GetView(row);
    return QString::fromUtf8(bytes.data(), int(bytes.size()));
}

template<typename Builder, typename Value>
arrow::Result<std::shared_ptr<arrow::Array>> array(int count, Value value) {
    Builder builder;
    ARROW_RETURN_NOT_OK(builder.Reserve(count));
    for (int row = 0; row < count; ++row) builder.UnsafeAppend(value(row));
    std::shared_ptr<arrow::Array> result;
    ARROW_RETURN_NOT_OK(builder.Finish(&result));
    return result;
}

arrow::Result<bool> writeCombined(
    const QVector<CandidateBundleIO::View> &views,
    const CandidatePoolRescorer::Result &result,
    const QVector<double> &probabilities, const QVector<int> &destinations,
    const QString &path, int familyFolds) {
    using namespace CandidateScoresReaderRowNamespace;
    auto fields = reportFields();
    std::sort(fields.begin(), fields.end(),
              [](const Field &a, const Field &b) { return a.name < b.name; });
    std::vector<std::shared_ptr<arrow::Table>> tables;
    QVector<int64_t> folds(result.confidence.inputIndices.size(), -1);
    int offset = 0;
    for (const auto &view : views) {
        const auto source = QDir(view.directory).filePath("report.parquet");
        if (view.reportSha256.size() != 64 || CandidateBundleIO::fileHash(source) != view.reportSha256)
            return arrow::Status::Invalid("Candidate report hash differs");
        ARROW_ASSIGN_OR_RAISE(auto table, readTable(source));
        if (table->num_rows() != view.candidates.identities.size())
            return arrow::Status::Invalid("Candidate report row count differs");
        if (!canonicalTable(*table, fields)) return false;
        const auto peptides = column<arrow::StringArray>(*table, PEP_STR_W_MODS);
        const auto origins = column<arrow::StringArray>(*table, PEP_STR_W_MODS_DECOY_OG);
        const auto targets = column<arrow::StringArray>(*table, TARG_KEY);
        const auto charges = column<arrow::FloatArray>(*table, CHARGE);
        const auto decoys = column<arrow::Int64Array>(*table, IS_DECOY);
        const auto apexes = column<arrow::FloatArray>(*table, SCAN_TIME);
        for (int row = 0; row < table->num_rows(); ++row) {
            const auto &identity = view.candidates.identities[row];
            if (text(*peptides, row) != identity.reportedPeptide
                || text(*origins, row) != identity.originPeptide
                || text(*targets, row) != identity.targetKey
                || charges->Value(row) != identity.charge
                || bool(decoys->Value(row)) != identity.isDecoy
                || apexes->Value(row) != identity.apex)
                return arrow::Status::Invalid("Candidate report identity differs");
            const int destination = destinations[offset + row];
            if (destination >= 0)
                folds[destination] = PeptideFamilyNeuralNet::familyHash(identity.originPeptide) % familyFolds;
        }
        offset += int(table->num_rows());
        tables.push_back(std::move(table));
    }
    const int count = result.confidence.inputIndices.size();
    if (count == 0) return arrow::Status::Invalid("No retained candidate rows");
    ARROW_ASSIGN_OR_RAISE(auto indices, array<arrow::Int64Builder>(count,
        [&](int row) { return int64_t(result.confidence.inputIndices[row]); }));
    ARROW_ASSIGN_OR_RAISE(auto all, arrow::ConcatenateTables(tables));
    ARROW_ASSIGN_OR_RAISE(auto selected, arrow::compute::Take(arrow::Datum(all), arrow::Datum(indices)));
    std::vector<std::shared_ptr<arrow::Field>> canonicalFields;
    for (const auto &field : fields)
        canonicalFields.push_back(arrow::field(field.name.toStdString(), field.type));
    // The row writer discards input schema metadata, so do the same here.
    auto output = arrow::Table::Make(arrow::schema(canonicalFields), selected.table()->columns());
    std::vector<std::pair<QString, std::shared_ptr<arrow::Array>>> replacements;
    ARROW_ASSIGN_OR_RAISE(auto scores, array<arrow::DoubleBuilder>(count,
        [&](int row) { return probabilities[result.confidence.inputIndices[row]]; }));
    replacements.push_back({CLASS_SCR, scores});
    ARROW_ASSIGN_OR_RAISE(auto foldArray, array<arrow::Int64Builder>(count,
        [&](int row) { return folds[row]; }));
    replacements.push_back({CLASS_FOLD, foldArray});
    ARROW_ASSIGN_OR_RAISE(auto qValues, array<arrow::DoubleBuilder>(count,
        [&](int row) { return result.confidence.qValues[row]; }));
    replacements.push_back({Q_VAL, qValues});
    replacements.push_back({PRECURSOR_Q_VAL, qValues});
    ARROW_ASSIGN_OR_RAISE(auto unassigned, array<arrow::DoubleBuilder>(count, [](int) { return 1.0; }));
    replacements.push_back({PEPTIDE_Q_VAL, unassigned});
    replacements.push_back({PROTEIN_Q_VAL, unassigned});
    ARROW_ASSIGN_OR_RAISE(auto noRatio, array<arrow::DoubleBuilder>(count, [](int) { return -1.0; }));
    replacements.push_back({DECOY_RATIO, noRatio});
    ARROW_ASSIGN_OR_RAISE(auto best, array<arrow::Int64Builder>(count, [](int) { return int64_t(1); }));
    replacements.push_back({IS_BEST_PRECURSOR_CANDIDATE, best});
    ARROW_ASSIGN_OR_RAISE(auto notBest, array<arrow::Int64Builder>(count, [](int) { return int64_t(0); }));
    replacements.push_back({IS_BEST_PEPTIDE_CANDIDATE, notBest});
    replacements.push_back({IS_BEST_PROTEIN_CANDIDATE, notBest});
    for (const auto &replacement : replacements) {
        const int index = output->schema()->GetFieldIndex(replacement.first.toStdString());
        ARROW_ASSIGN_OR_RAISE(output, output->SetColumn(
            index, output->field(index), std::make_shared<arrow::ChunkedArray>(replacement.second)));
    }
    ARROW_RETURN_NOT_OK(writeTable(*output, path));
    return true;
}
} // namespace

Error::Err CandidateReportColumns::write(
    const QVector<CandidateScoresReaderRow> &rows, const QString &path) {
    if (rows.isEmpty() || path.isEmpty() || QFile::exists(path)) return Error::eValueError;
    const auto status = writeRows(rows, path);
    if (!status.ok()) {
        qWarning() << "Could not write candidate report:" << status.ToString().c_str();
        return Error::eError;
    }
    return Error::eNoError;
}

Error::Err CandidateReportColumns::tryWriteCombined(
    const QVector<CandidateBundleIO::View> &views,
    const CandidatePoolRescorer::Result &result,
    const QVector<double> &probabilities, const QVector<int> &destinations,
    const QString &path, int familyFolds, bool *written) {
    if (written == nullptr) return Error::eValueError;
    *written = false;
    const auto outcome = writeCombined(views, result, probabilities, destinations, path, familyFolds);
    if (!outcome.ok()) {
        qWarning() << "Could not combine candidate report columns:" << outcome.status().ToString().c_str();
        return Error::eValueError;
    }
    *written = *outcome;
    return Error::eNoError;
}
