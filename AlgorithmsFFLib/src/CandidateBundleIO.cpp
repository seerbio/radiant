#include "CandidateBundleIO.h"
#include "CandidateFeatureSchema.h"
#include "CandidateReportColumns.h"
#include "ParquetReader.h"

#include <QCryptographicHash>
#include <QDataStream>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QSet>
#include <QSysInfo>
#include <cmath>
#include <cstring>

namespace {
constexpr quint32 maximumRows = 1000000;
const QString metadataName = QStringLiteral("candidates.bin");
const QString featuresName = QStringLiteral("features.f32");
const QString reportName = QStringLiteral("report.parquet");
bool newDirectory(const QString &path) {
    const QFileInfo info(path);
    QDir parent(info.absolutePath());
    return QDir().mkpath(parent.absolutePath()) && parent.mkdir(info.fileName());
}
void writeString(QDataStream &stream, const QString &text) {
    const auto bytes = text.toUtf8();
    stream << quint32(bytes.size());
    stream.writeRawData(bytes.constData(), bytes.size());
}
bool readString(QDataStream &stream, QString *text) {
    quint32 size = 0;
    stream >> size;
    if (stream.status() != QDataStream::Ok || size > 1048576) return false;
    QByteArray bytes(int(size), Qt::Uninitialized);
    if (stream.readRawData(bytes.data(), int(size)) != int(size)) return false;
    *text = QString::fromUtf8(bytes);
    return text->toUtf8() == bytes;
}
FragmentCompetition::Evidence evidence(const CandidateScoresReaderRow &row) {
    FragmentCompetition::Evidence result;
    result.reportedPeptide = row.peptideStringWithMods;
    result.mass = row.mass;
    result.charge = row.charge;
    result.apex = row.scanTime;
    result.width = double(row.scanTimeEnd) - row.scanTimeStart;
    result.isDecoy = row.isDecoy;
    result.priority = row.discriminantScore;
    result.searched = {row.mzSearched1, row.mzSearched2, row.mzSearched3, row.mzSearched4,
                       row.mzSearched5, row.mzSearched6, row.mzSearched7, row.mzSearched8,
                       row.mzSearched9, row.mzSearched10, row.mzSearched11, row.mzSearched12};
    result.intensity = {row.intensityFoundMax1, row.intensityFoundMax2, row.intensityFoundMax3,
                        row.intensityFoundMax4, row.intensityFoundMax5, row.intensityFoundMax6,
                        row.intensityFoundMax7, row.intensityFoundMax8, row.intensityFoundMax9,
                        row.intensityFoundMax10, row.intensityFoundMax11, row.intensityFoundMax12};
    result.cosine = {row.cosineSimToAnchor1, row.cosineSimToAnchor2, row.cosineSimToAnchor3,
                     row.cosineSimToAnchor4, row.cosineSimToAnchor5, row.cosineSimToAnchor6,
                     row.cosineSimToAnchor7, row.cosineSimToAnchor8, row.cosineSimToAnchor9,
                     row.cosineSimToAnchor10, row.cosineSimToAnchor11, row.cosineSimToAnchor12};
    return result;
}
bool validIdentity(const CandidatePoolSelection::Identity &id, const FragmentCompetition::Evidence &e) {
    return !id.reportedPeptide.isEmpty() && !id.originPeptide.isEmpty() && !id.targetKey.isEmpty()
        && id.charge > 0 && std::isfinite(e.mass) && e.mass > 0
        && std::isfinite(e.charge) && id.charge == e.charge && std::isfinite(e.apex)
        && std::isfinite(e.width) && e.width > 0;
}
}

QString CandidateBundleIO::fileHash(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return {};
    QCryptographicHash hash(QCryptographicHash::Sha256);
    if (!hash.addData(&file) || file.error() != QFileDevice::NoError) return {};
    return QString::fromLatin1(hash.result().toHex());
}

Error::Err CandidateBundleIO::write(
    const QVector<CandidateScores*> &rows, const QString &runName,
    const QJsonObject &provenance, const QString &directory) {
    if (QSysInfo::ByteOrder != QSysInfo::LittleEndian || rows.isEmpty()
        || rows.size() > maximumRows || runName.isEmpty() || directory.isEmpty()
        || QFile::exists(directory)) return Error::eValueError;
    for (const auto *row : rows) {
        if (row == nullptr || row->targetDecoyCandidatePair == nullptr
            || row->featuresArray.size() != FeaturesSize || row->integrations.size() < 12
            || row->ionLabels.size() < 12 || !std::isfinite(row->discriminantScore))
            return Error::eValueError;
    }
    if (!newDirectory(directory)) return Error::eFileError;
    const QDir dir(directory);
    QFile metadata(dir.filePath(metadataName)), features(dir.filePath(featuresName));
    if (!metadata.open(QIODevice::WriteOnly | QIODevice::NewOnly)
        || !features.open(QIODevice::WriteOnly | QIODevice::NewOnly)) return Error::eFileError;
    if (metadata.write("RDCB0001", 8) != 8) return Error::eFileError;
    QDataStream stream(&metadata);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    stream << quint32(FeaturesSize) << quint32(rows.size());
    QVector<CandidateScoresReaderRow> report;
    report.reserve(rows.size());
    for (int index = 0; index < rows.size(); ++index) {
        const auto *source = rows[index];
        auto row = CandidateScoresReaderRow::buildCandidateScoresReaderRow(source);
        const auto e = evidence(row);
        if (!std::isfinite(row.charge) || row.charge <= 0 || row.charge > INT_MAX
            || row.charge != std::floor(row.charge)) return Error::eValueError;
        CandidatePoolSelection::Identity id{
            row.peptideStringWithMods, row.peptideStringWithModsDecoyOrigin, row.targetKey,
            int(row.charge), 0, index, row.isDecoy, row.scanTime};
        if (!validIdentity(id, e)) return Error::eValueError;
        writeString(stream, id.reportedPeptide);
        writeString(stream, id.originPeptide);
        writeString(stream, id.targetKey);
        stream << e.charge << e.mass << e.apex << e.width << quint32(e.isDecoy);
        for (double value : e.searched) stream << value;
        for (double value : e.intensity) stream << value;
        const qint64 size = FeaturesSize * sizeof(float);
        if (features.write(reinterpret_cast<const char*>(source->featuresArray.constData()), size) != size)
            return Error::eFileError;
        report.push_back(std::move(row));
    }
    if (stream.status() != QDataStream::Ok || !metadata.flush() || !features.flush())
        return Error::eFileError;
    metadata.close();
    features.close();
    const auto error = CandidateReportColumns::write(report, dir.filePath(reportName));
    if (error != Error::eNoError) return error;
    QJsonObject files;
    for (const auto &name : {metadataName, featuresName, reportName}) {
        const auto hash = fileHash(dir.filePath(name));
        if (hash.isEmpty()) return Error::eFileError;
        files[name] = hash;
    }
    const QJsonObject manifest{
        {"format", "radiant-candidate-bundle"}, {"version", 1}, {"rows", rows.size()},
        {"feature_schema_sha256", CandidateFeatureSchema::id()}, {"feature_names", CandidateFeatureSchema::names()},
        {"run_name", runName}, {"files", files}, {"provenance", provenance}};
    QFile completion(dir.filePath("manifest.json"));
    if (!completion.open(QIODevice::WriteOnly | QIODevice::NewOnly)) return Error::eFileError;
    const auto bytes = QJsonDocument(manifest).toJson(QJsonDocument::Indented);
    if (completion.write(bytes) != bytes.size() || !completion.flush()) return Error::eFileError;
    return Error::eNoError;
}

Error::Err CandidateBundleIO::read(const QString &directory, View *view) {
    if (view == nullptr || QSysInfo::ByteOrder != QSysInfo::LittleEndian) return Error::eValueError;
    const QDir dir(directory);
    QFile manifestFile(dir.filePath("manifest.json"));
    if (!manifestFile.open(QIODevice::ReadOnly)) return Error::eFileError;
    QJsonParseError parseError;
    const auto manifestBytes = manifestFile.readAll();
    const auto document = QJsonDocument::fromJson(manifestBytes, &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) return Error::eValueError;
    const auto manifest = document.object();
    const QSet<QString> allowed{"format", "version", "rows", "feature_schema_sha256", "feature_names",
                                "run_name", "files", "provenance"};
    for (const auto &key : manifest.keys()) if (!allowed.contains(key)) return Error::eValueError;
    const int count = manifest["rows"].toInt(-1);
    if (manifest["format"].toString() != "radiant-candidate-bundle" || manifest["version"].toInt(-1) != 1
        || count < 1 || count > maximumRows || manifest["run_name"].toString().isEmpty()
        || manifest["feature_schema_sha256"].toString() != CandidateFeatureSchema::id()
        || manifest["feature_names"].toArray() != CandidateFeatureSchema::names()
        || !manifest["provenance"].isObject() || !manifest["files"].isObject()) return Error::eValueError;
    const auto hashes = manifest["files"].toObject();
    if (hashes.size() != 3) return Error::eValueError;
    for (const auto &name : {metadataName, featuresName, reportName}) {
        const QString expected = hashes[name].toString();
        if (expected.size() != 64 || fileHash(dir.filePath(name)) != expected) return Error::eValueError;
    }
    QFile metadata(dir.filePath(metadataName)), featureFile(dir.filePath(featuresName));
    if (!metadata.open(QIODevice::ReadOnly) || metadata.read(8) != "RDCB0001"
        || !featureFile.open(QIODevice::ReadOnly)
        || featureFile.size() != qint64(count) * FeaturesSize * sizeof(float)) return Error::eValueError;
    QDataStream stream(&metadata);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    quint32 width = 0, rows = 0;
    stream >> width >> rows;
    if (width != FeaturesSize || rows != quint32(count)) return Error::eValueError;
    View result;
    result.directory = dir.absolutePath();
    result.reportSha256 = hashes[reportName].toString();
    result.manifestSha256 = QString::fromLatin1(
        QCryptographicHash::hash(manifestBytes, QCryptographicHash::Sha256).toHex());
    result.fileSha256 = hashes;
    result.provenance = manifest["provenance"].toObject();
    auto &run = result.candidates;
    run.name = manifest["run_name"].toString();
    run.identities.reserve(count);
    run.evidence.reserve(count);
    run.rawFeatures.reserve(count);
    for (int index = 0; index < count; ++index) {
        CandidatePoolSelection::Identity id;
        FragmentCompetition::Evidence e;
        quint32 label = 0;
        if (!readString(stream, &id.reportedPeptide) || !readString(stream, &id.originPeptide)
            || !readString(stream, &id.targetKey)) return Error::eValueError;
        stream >> e.charge >> e.mass >> e.apex >> e.width >> label;
        if (label > 1 || !std::isfinite(e.charge) || e.charge <= 0 || e.charge > INT_MAX
            || std::floor(e.charge) != e.charge) return Error::eValueError;
        id.charge = int(e.charge);
        id.isDecoy = e.isDecoy = label;
        id.apex = e.apex;
        id.viewRank = index;
        e.reportedPeptide = id.reportedPeptide;
        for (double &value : e.searched) stream >> value;
        for (double &value : e.intensity) stream >> value;
        QVector<float> features(FeaturesSize);
        const qint64 size = FeaturesSize * sizeof(float);
        if (featureFile.read(reinterpret_cast<char*>(features.data()), size) != size
            || stream.status() != QDataStream::Ok || !validIdentity(id, e)
            || features[Mass] != e.mass || features[Charge] != e.charge) return Error::eValueError;
        for (int ion = 0; ion < 12; ++ion) e.cosine[ion] = features[CosineSimToAnchor1 + ion];
        run.identities.push_back(std::move(id));
        run.evidence.push_back(std::move(e));
        run.rawFeatures.push_back(std::move(features));
    }
    if (!metadata.atEnd() || !featureFile.atEnd()) return Error::eValueError;
    *view = std::move(result);
    return Error::eNoError;
}

Error::Err CandidateBundleIO::writeReport(
    const QVector<View> &views, const CandidatePoolRescorer::Result &result, const QString &path, int familyFolds) {
    if (views.isEmpty() || path.isEmpty() || QFile::exists(path) || familyFolds < 2
        || result.selectedInputIndices.size() != result.combinedProbability.size()
        || result.confidence.inputIndices.size() != result.confidence.qValues.size())
        return Error::eValueError;
    qint64 total = 0;
    for (const auto &view : views) total += view.candidates.identities.size();
    if (total > INT_MAX) return Error::eValueError;
    QVector<double> probabilities(int(total), -1);
    QSet<int> seen;
    for (int row = 0; row < result.selectedInputIndices.size(); ++row) {
        const int index = result.selectedInputIndices[row];
        const double p = result.combinedProbability[row];
        if (index < 0 || index >= total || seen.contains(index)
            || !std::isfinite(p) || p < 0 || p > 1) return Error::eValueError;
        seen.insert(index);
        probabilities[index] = p;
    }
    QVector<int> destinations(int(total), -1);
    for (int row = 0; row < result.confidence.inputIndices.size(); ++row) {
        const int index = result.confidence.inputIndices[row];
        const double q = result.confidence.qValues[row];
        if (index < 0 || index >= total || destinations[index] != -1 || probabilities[index] < 0
            || !std::isfinite(q) || q < 0 || q > 1) return Error::eValueError;
        destinations[index] = row;
    }
    bool written = false;
    const auto columnError = CandidateReportColumns::tryWriteCombined(
        views, result, probabilities, destinations, path, familyFolds, &written);
    if (columnError != Error::eNoError || written) return columnError;
    // Older or noncanonical reports retain the original typed row conversion.
    QVector<CandidateScoresReaderRow> output(result.confidence.inputIndices.size());
    int offset = 0;
    for (const auto &view : views) {
        // Bind the report to the manifest verified when the view was read.
        const QString reportPath = QDir(view.directory).filePath(reportName);
        if (view.reportSha256.size() != 64 || view.reportSha256 != fileHash(reportPath))
            return Error::eValueError;
        QVector<CandidateScoresReaderRow> report;
        const auto error = ParquetReader::read(reportPath, &report);
        if (error != Error::eNoError) return error;
        if (report.size() != view.candidates.identities.size()) return Error::eValueError;
        for (int row = 0; row < report.size(); ++row) {
            const auto &id = view.candidates.identities[row];
            auto &record = report[row];
            if (record.peptideStringWithMods != id.reportedPeptide
                || record.peptideStringWithModsDecoyOrigin != id.originPeptide || record.targetKey != id.targetKey
                || record.charge != id.charge || record.isDecoy != id.isDecoy || record.scanTime != id.apex)
                return Error::eValueError;
            const int destination = destinations[offset + row];
            if (destination < 0) continue;
            record.classifierScore = probabilities[offset + row];
            record.decoyRatio = -1; // The exported LDA diagnostic is not a combined-model confidence estimate.
            record.classifierFold = int(PeptideFamilyNeuralNet::familyHash(id.originPeptide) % familyFolds);
            record.qValue = record.precursorQValue = result.confidence.qValues[destination];
            record.isBestPrecursorCandidate = 1;
            record.peptideQValue = record.proteinQValue = 1;
            record.isBestPeptideCandidate = record.isBestProteinCandidate = 0;
            output[destination] = std::move(record);
        }
        offset += report.size();
    }
    if (output.isEmpty()) return Error::eValueError;
    return CandidateReportColumns::write(output, path);
}
