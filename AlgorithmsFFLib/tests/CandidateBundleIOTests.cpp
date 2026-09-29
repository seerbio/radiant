#include "CandidateBundleIO.h"
#include "CandidateFeatureSchema.h"
#include "CandidateReportColumns.h"
#include "FragmentCompetitionTestData.h"

#include <QFile>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <QtTest>
#include <cstring>

class CandidateBundleIOTests : public QObject {
    Q_OBJECT
private slots:
    void everyReportFieldRoundTrips();
    void reconstructedFeaturesUseMatchingColumns();
    void roundTripAndReport();
    void rejectsCorruption();
    void invalidInputDoesNotOverwrite();
    void columnWriterMatchesLegacyBits();
    void combinedColumnsMatchLegacyAcrossViews();
    void noncanonicalReportRetainsLegacyConversion_data();
    void noncanonicalReportRetainsLegacyConversion();
    void rejectsIdentityMismatchInUnselectedRow();
};

namespace {
void prepare(CompetitionCandidateFixture &fixture, bool decoy) {
    fixture.scores.isDecoy = decoy;
    fixture.scores.targetKey = "500123";
    fixture.scores.proteinGroup = decoy ? "decoy_annotation" : "target_annotation";
    fixture.scores.featuresArray[CosineSim45MS1] = .45f;
    fixture.scores.featuresArray[CosineSimSum45] = .79f;
    fixture.scores.featuresArray[IntensityFoundMax1] = .125f;
    fixture.scores.classifierScore = .8;
    fixture.scores.decoyRatio = .37;
    fixture.scores.qValue = .9;
    fixture.scores.precursorQValue = .9;
    fixture.scores.peptideQValue = .01;
    fixture.scores.proteinQValue = .02;
    fixture.scores.isBestPeptideCandidate = 1;
    fixture.scores.isBestProteinCandidate = 1;
}
QByteArray readBytes(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return {};
    return file.readAll();
}
bool replaceBytes(const QString &path, const QByteArray &bytes) {
    QFile file(path);
    return file.open(QIODevice::WriteOnly | QIODevice::Truncate)
        && file.write(bytes) == bytes.size() && file.flush();
}
struct MappedRow : ParquetReaderInputBase {
    QMap<QString, QVariant> map() override { return dataMap(); }
};
bool updateReportHash(const QString &directory) {
    const QString path = directory + "/manifest.json";
    auto manifest = QJsonDocument::fromJson(readBytes(path)).object();
    auto files = manifest["files"].toObject();
    files["report.parquet"] = CandidateBundleIO::fileHash(directory + "/report.parquet");
    manifest["files"] = files;
    return replaceBytes(path, QJsonDocument(manifest).toJson());
}
QVector<CandidateScoresReaderRow> legacyCombined(
    const QVector<CandidateBundleIO::View> &views, const CandidatePoolRescorer::Result &result) {
    QVector<CandidateScoresReaderRow> all;
    for (const auto &view : views) {
        QVector<CandidateScoresReaderRow> rows;
        if (ParquetReader::read(view.directory + "/report.parquet", &rows) != Error::eNoError)
            return {};
        all += rows;
    }
    QVector<CandidateScoresReaderRow> output;
    for (int row = 0; row < result.confidence.inputIndices.size(); ++row) {
        const int index = result.confidence.inputIndices[row];
        auto record = all[index];
        record.classifierScore = result.combinedProbability[result.selectedInputIndices.indexOf(index)];
        record.decoyRatio = -1;
        record.classifierFold = int(PeptideFamilyNeuralNet::familyHash(
            record.peptideStringWithModsDecoyOrigin) % 3);
        record.qValue = record.precursorQValue = result.confidence.qValues[row];
        record.isBestPrecursorCandidate = 1;
        record.peptideQValue = record.proteinQValue = 1;
        record.isBestPeptideCandidate = record.isBestProteinCandidate = 0;
        output.push_back(std::move(record));
    }
    return output;
}
}

void CandidateBundleIOTests::columnWriterMatchesLegacyBits() {
    CandidateScoresReaderRow defaults;
    QVector<CandidateScoresReaderRow> rows;
    for (int row = 0; row < 37; ++row) {
        auto values = defaults.map();
        int ordinal = 0;
        for (auto it = values.begin(); it != values.end(); ++it) {
            ++ordinal;
            switch (it.value().userType()) {
            case QMetaType::Float: it.value() = float(ordinal) + float(row) * .125f; break;
            case QMetaType::Double: it.value() = double(ordinal) + double(row) * .125; break;
            case QMetaType::Int: it.value() = ordinal * (row % 2 ? -1 : 1); break;
            case QMetaType::Bool: it.value() = bool(row % 2); break;
            case QMetaType::QString:
                it.value() = QString::fromUtf8("λ_") + it.key() + QString::number(row); break;
            default: QFAIL("Unhandled field");
            }
        }
        values["PeptideStringWithMods"] = "_PEPTIDEK_";
        values["PeptideStringWithModsDecoyOrigin"] = "_ASDFGHK_";
        ParquetReaderInputBase mapped;
        mapped.setDataMap(values);
        CandidateScoresReaderRow record;
        QCOMPARE(record.initFromRead(mapped), Error::eNoError);
        rows.push_back(std::move(record));
    }
    rows[0].mass = -0.0f;
    rows[0].classifierScore = -0.0;
    const quint32 floatNaN = 0x7fc01234;
    const quint64 doubleNaN = 0x7ff8000000001234ULL;
    std::memcpy(&rows[1].mass, &floatNaN, sizeof(float));
    std::memcpy(&rows[1].classifierScore, &doubleNaN, sizeof(double));
    QTemporaryDir directory;
    const auto legacy = directory.filePath("legacy.parquet");
    const auto columnar = directory.filePath("columns.parquet");
    QCOMPARE(ParquetReader::write(rows, legacy), Error::eNoError);
    QCOMPARE(CandidateReportColumns::write(rows, columnar), Error::eNoError);
    QVERIFY(!readBytes(legacy).isEmpty());
    QCOMPARE(readBytes(columnar), readBytes(legacy));
    QVERIFY(CandidateReportColumns::write(rows, columnar) != Error::eNoError);
}

void CandidateBundleIOTests::combinedColumnsMatchLegacyAcrossViews() {
    QTemporaryDir directory;
    CompetitionCandidateFixture first("PEPTIDEK"), second("AGVTFERK");
    CompetitionCandidateFixture third("PEPTIDER"), fourth("GILGFVFTLK");
    prepare(first, false);
    prepare(second, true);
    prepare(third, false);
    prepare(fourth, true);
    const auto one = directory.filePath("one");
    const auto two = directory.filePath("two");
    QCOMPARE(CandidateBundleIO::write({&first.scores, &second.scores}, "sample", {}, one), Error::eNoError);
    QCOMPARE(CandidateBundleIO::write({&third.scores, &fourth.scores}, "sample", {}, two), Error::eNoError);
    CandidateBundleIO::View left, right;
    QCOMPARE(CandidateBundleIO::read(one, &left), Error::eNoError);
    QCOMPARE(CandidateBundleIO::read(two, &right), Error::eNoError);
    CandidatePoolRescorer::Result result;
    result.selectedInputIndices = {3, 0, 2, 1};
    result.combinedProbability = {.1, .3, .2, .7};
    result.confidence.inputIndices = {2, 0, 3};
    result.confidence.qValues = {.01, .02, .4};
    const auto expected = legacyCombined({left, right}, result);
    QCOMPARE(expected.size(), 3);
    const auto legacy = directory.filePath("legacy.parquet");
    const auto columnar = directory.filePath("columns.parquet");
    QCOMPARE(ParquetReader::write(expected, legacy), Error::eNoError);
    QCOMPARE(CandidateBundleIO::writeReport({left, right}, result, columnar), Error::eNoError);
    QCOMPARE(readBytes(columnar), readBytes(legacy));
}

void CandidateBundleIOTests::noncanonicalReportRetainsLegacyConversion_data() {
    QTest::addColumn<bool>("extraColumn");
    QTest::newRow("extra-column") << true;
    QTest::newRow("noncanonical-boolean") << false;
}

void CandidateBundleIOTests::noncanonicalReportRetainsLegacyConversion() {
    QFETCH(bool, extraColumn);
    QTemporaryDir directory;
    CompetitionCandidateFixture candidate;
    prepare(candidate, true);
    const auto viewPath = directory.filePath("view");
    QCOMPARE(CandidateBundleIO::write({&candidate.scores}, "sample", {}, viewPath), Error::eNoError);
    auto values = CandidateScoresReaderRow::buildCandidateScoresReaderRow(&candidate.scores).map();
    if (extraColumn) values["UnusedAnnotation"] = "ignored by the original row conversion";
    else values["IsDecoy"] = 7;
    MappedRow row;
    row.setDataMap(values);
    QCOMPARE(ParquetReader::write(QVector<MappedRow>{row}, viewPath + "/report.parquet"), Error::eNoError);
    QVERIFY(updateReportHash(viewPath));
    CandidateBundleIO::View view;
    QCOMPARE(CandidateBundleIO::read(viewPath, &view), Error::eNoError);
    CandidatePoolRescorer::Result result;
    result.selectedInputIndices = {0};
    result.combinedProbability = {.125};
    result.confidence.inputIndices = {0};
    result.confidence.qValues = {.5};
    const auto expected = legacyCombined({view}, result);
    QCOMPARE(expected.size(), 1);
    const auto legacy = directory.filePath("legacy.parquet");
    const auto output = directory.filePath("output.parquet");
    QCOMPARE(ParquetReader::write(expected, legacy), Error::eNoError);
    QCOMPARE(CandidateBundleIO::writeReport({view}, result, output), Error::eNoError);
    QCOMPARE(readBytes(output), readBytes(legacy));
}

void CandidateBundleIOTests::rejectsIdentityMismatchInUnselectedRow() {
    QTemporaryDir directory;
    CompetitionCandidateFixture first("PEPTIDEK"), second("ANOTHERK");
    prepare(first, false);
    prepare(second, true);
    const auto viewPath = directory.filePath("view");
    QCOMPARE(CandidateBundleIO::write({&first.scores, &second.scores}, "sample", {}, viewPath), Error::eNoError);
    QVector<CandidateScoresReaderRow> rows;
    QCOMPARE(ParquetReader::read(viewPath + "/report.parquet", &rows), Error::eNoError);
    rows[1].targetKey = "different";
    QCOMPARE(ParquetReader::write(rows, viewPath + "/report.parquet"), Error::eNoError);
    QVERIFY(updateReportHash(viewPath));
    CandidateBundleIO::View view;
    QCOMPARE(CandidateBundleIO::read(viewPath, &view), Error::eNoError);
    CandidatePoolRescorer::Result result;
    result.selectedInputIndices = {0};
    result.combinedProbability = {.125};
    result.confidence.inputIndices = {0};
    result.confidence.qValues = {.5};
    const auto output = directory.filePath("must-not-exist.parquet");
    QVERIFY(CandidateBundleIO::writeReport({view}, result, output) != Error::eNoError);
    QVERIFY(!QFile::exists(output));
}

void CandidateBundleIOTests::everyReportFieldRoundTrips() {
    CandidateScoresReaderRow defaults;
    auto expected = defaults.map();
    int ordinal = 0;
    for (auto it = expected.begin(); it != expected.end(); ++it) {
        ++ordinal;
        // Distinct exactly representable values expose omitted and crossed
        // columns, including columns that normal fixtures leave at zero.
        switch (it.value().userType()) {
        case QMetaType::Float: it.value() = float(ordinal) + .125f; break;
        case QMetaType::Double: it.value() = double(ordinal) + .125; break;
        case QMetaType::Int: it.value() = ordinal; break;
        case QMetaType::Bool: it.value() = true; break;
        case QMetaType::QString: it.value() = "value_" + it.key(); break;
        default: QFAIL(qPrintable("Unhandled report type: " + it.key()));
        }
    }
    expected["PeptideStringWithMods"] = "_PEPTIDEK_";
    expected["PeptideStringWithModsDecoyOrigin"] = "_ASDFGHK_";
    ParquetReaderInputBase input;
    input.setDataMap(expected);
    CandidateScoresReaderRow decoded;
    QCOMPARE(decoded.initFromRead(input), Error::eNoError);
    QCOMPARE(decoded.map(), expected);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath("all-fields.parquet");
    QCOMPARE(ParquetReader::write(QVector<CandidateScoresReaderRow>{decoded}, path), Error::eNoError);
    QVector<CandidateScoresReaderRow> reloaded;
    QCOMPARE(ParquetReader::read(path, &reloaded), Error::eNoError);
    QCOMPARE(reloaded.size(), 1);
    QCOMPARE(reloaded[0].map(), expected);
}

void CandidateBundleIOTests::reconstructedFeaturesUseMatchingColumns() {
    CandidateScoresReaderRow row;
    row.cosineSim45MS1 = .45f;
    row.cosineSimSum45 = .79f;
    row.altTargetKeyIdDiscScoreCharge3_1 = .31f;
    row.altTargetKeyIdDiscScoreCharge3_2 = .32f;
    row.altTargetKeyIdDiscScoreCharge4_1 = .41f;
    row.altTargetKeyIdDiscScoreCharge4_2 = .42f;
    const auto features = CandidateScoresReaderRow::featuresArrayFromCandidateScoresReaderRow(row);
    QCOMPARE(features[CosineSim45MS1], .45f);
    QCOMPARE(features[CosineSimSum45], .79f);
    QCOMPARE(features[AltTargetKeyIdDiscScoreCharge3_1], .31f);
    QCOMPARE(features[AltTargetKeyIdDiscScoreCharge3_2], .32f);
    QCOMPARE(features[AltTargetKeyIdDiscScoreCharge4_1], .41f);
    QCOMPARE(features[AltTargetKeyIdDiscScoreCharge4_2], .42f);
}

void CandidateBundleIOTests::roundTripAndReport() {
    QTemporaryDir temp;
    QVERIFY(temp.isValid());
    CompetitionCandidateFixture first("PEPTIDEK"), second("ANOTHERK");
    prepare(first, false);
    prepare(second, true);
    QVector<CandidateScores*> candidates{&first.scores, &second.scores};
    const QString path = temp.filePath("view");
    const QJsonObject provenance{{"source_sha256", QString(64, 'a')}, {"library_sha256", QString(64, 'b')}};
    QCOMPARE(CandidateBundleIO::write(candidates, "sample", provenance, path), Error::eNoError);
    CandidateBundleIO::View loaded;
    QCOMPARE(CandidateBundleIO::read(path, &loaded), Error::eNoError);
    QCOMPARE(loaded.candidates.name, QString("sample"));
    QCOMPARE(loaded.provenance, provenance);
    QCOMPARE(loaded.manifestSha256, CandidateBundleIO::fileHash(path + "/manifest.json"));
    QCOMPARE(loaded.fileSha256["report.parquet"].toString(), loaded.reportSha256);
    QCOMPARE(loaded.candidates.rawFeatures.size(), 2);
    for (int i = 0; i < candidates.size(); ++i) {
        QVERIFY(std::memcmp(candidates[i]->featuresArray.constData(),
                            loaded.candidates.rawFeatures[i].constData(), FeaturesSize * sizeof(float)) == 0);
        const auto report = CandidateScoresReaderRow::buildCandidateScoresReaderRow(candidates[i]);
        const auto &id = loaded.candidates.identities[i];
        const auto &e = loaded.candidates.evidence[i];
        QCOMPARE(id.reportedPeptide, QString(report.peptideStringWithMods));
        QCOMPARE(id.originPeptide, QString(report.peptideStringWithModsDecoyOrigin));
        QCOMPARE(id.isDecoy, report.isDecoy);
        QCOMPARE(id.viewRank, i);
        QCOMPARE(e.width, double(report.scanTimeEnd) - report.scanTimeStart);
        QCOMPARE(e.intensity[0], double(candidates[i]->integrations[0]));
        QCOMPARE(e.searched[0], double(report.mzSearched1));
        QCOMPARE(e.cosine[0], double(candidates[i]->featuresArray[CosineSimToAnchor1]));
        QCOMPARE(loaded.candidates.rawFeatures[i][IntensityFoundMax1], .125f);
        QCOMPARE(loaded.candidates.rawFeatures[i][CosineSim45MS1], .45f);
        QCOMPARE(loaded.candidates.rawFeatures[i][CosineSimSum45], .79f);
    }
    CandidatePoolRescorer::Result result;
    result.name = "sample";
    result.selectedInputIndices = {0, 1};
    result.combinedProbability = {.2, .1};
    result.confidence.inputIndices = {1, 0};
    result.confidence.qValues = {.6, .6};
    const QString output = temp.filePath("final.radiantDIA");
    QCOMPARE(CandidateBundleIO::writeReport({loaded}, result, output), Error::eNoError);
    QVector<CandidateScoresReaderRow> records;
    QCOMPARE(ParquetReader::read(output, &records), Error::eNoError);
    QCOMPARE(records.size(), 2);
    const QSet<QString> changed{
        "ClassifierScore", "ClassifierFold", "QValue", "PrecursorQValue", "PeptideQValue", "ProteinQValue",
        "IsBestPrecursorCandidate", "IsBestPeptideCandidate", "IsBestProteinCandidate", "DecoyRatio"};
    for (int i = 0; i < records.size(); ++i) {
        auto expected = CandidateScoresReaderRow::buildCandidateScoresReaderRow(candidates[1 - i]).map();
        auto actual = records[i].map();
        for (const auto &key : changed) { expected.remove(key); actual.remove(key); }
        for (const auto &key : expected.keys())
            if (actual.value(key) != expected.value(key))
                qWarning() << "Report round-trip mismatch" << key << actual.value(key) << expected.value(key);
        QCOMPARE(actual, expected);
        QCOMPARE(records[i].classifierScore, i == 0 ? .1 : .2);
        QCOMPARE(records[i].decoyRatio, -1.0);
        QCOMPARE(records[i].qValue, .6);
        QCOMPARE(records[i].precursorQValue, .6);
        QCOMPARE(records[i].isBestPrecursorCandidate, 1);
        QCOMPARE(records[i].peptideQValue, 1.0);
        QCOMPARE(records[i].proteinQValue, 1.0);
        QCOMPARE(records[i].isBestPeptideCandidate, 0);
        QCOMPARE(records[i].isBestProteinCandidate, 0);
    }
    // The original candidate scores and confidence remain untouched.
    QCOMPARE(first.scores.classifierScore, .8);
    QCOMPARE(first.scores.qValue, .9);
    QCOMPARE(first.scores.decoyRatio, .37);
}

void CandidateBundleIOTests::rejectsCorruption() {
    QTemporaryDir temp;
    CompetitionCandidateFixture first;
    prepare(first, false);
    const QString path = temp.filePath("view");
    QCOMPARE(CandidateBundleIO::write({&first.scores}, "sample", {}, path), Error::eNoError);
    CandidateBundleIO::View loaded;
    QCOMPARE(CandidateBundleIO::read(path, &loaded), Error::eNoError);
    const auto originalFeatures = readBytes(path + "/features.f32");
    auto corrupt = originalFeatures;
    corrupt[0] = char(corrupt[0] ^ 1);
    QVERIFY(replaceBytes(path + "/features.f32", corrupt));
    CandidateBundleIO::View untouched;
    untouched.candidates.name = "sentinel";
    QVERIFY(CandidateBundleIO::read(path, &untouched) != Error::eNoError);
    QCOMPARE(untouched.candidates.name, QString("sentinel"));
    QVERIFY(replaceBytes(path + "/features.f32", originalFeatures));
    const auto originalManifest = readBytes(path + "/manifest.json");
    const auto document = QJsonDocument::fromJson(originalManifest);
    for (const auto &field : {"version", "rows", "feature_schema_sha256", "feature_names"}) {
        auto manifest = document.object();
        if (QString(field) == "feature_schema_sha256") manifest[field] = QString(64, '0');
        else if (QString(field) == "feature_names") manifest[field] = QJsonArray{};
        else manifest[field] = 2;
        QVERIFY(replaceBytes(path + "/manifest.json", QJsonDocument(manifest).toJson()));
        QVERIFY(CandidateBundleIO::read(path, &untouched) != Error::eNoError);
        QCOMPARE(untouched.candidates.name, QString("sentinel"));
    }
    QVERIFY(replaceBytes(path + "/manifest.json", originalManifest));
    const auto originalReport = readBytes(path + "/report.parquet");
    QVERIFY(replaceBytes(path + "/report.parquet", originalReport + "bad"));
    CandidatePoolRescorer::Result result;
    result.selectedInputIndices = {0};
    result.combinedProbability = {.1};
    result.confidence.inputIndices = {0};
    result.confidence.qValues = {1};
    const auto output = temp.filePath("must-not-exist.radiantDIA");
    QVERIFY(CandidateBundleIO::writeReport({loaded}, result, output) != Error::eNoError);
    QVERIFY(!QFile::exists(output));
}

void CandidateBundleIOTests::invalidInputDoesNotOverwrite() {
    QTemporaryDir temp;
    CompetitionCandidateFixture first;
    prepare(first, false);
    const QString path = temp.filePath("view");
    QCOMPARE(CandidateBundleIO::write({&first.scores}, "sample", {}, path), Error::eNoError);
    const auto original = CandidateBundleIO::fileHash(path + "/manifest.json");
    QVERIFY(CandidateBundleIO::write({&first.scores}, "sample", {}, path) != Error::eNoError);
    QCOMPARE(CandidateBundleIO::fileHash(path + "/manifest.json"), original);
    const QString invalid = temp.filePath("invalid");
    QVERIFY(CandidateBundleIO::write({nullptr}, "sample", {}, invalid) != Error::eNoError);
    QVERIFY(!QFile::exists(invalid));
    CandidateBundleIO::View loaded;
    QCOMPARE(CandidateBundleIO::read(path, &loaded), Error::eNoError);
    CandidatePoolRescorer::Result result;
    result.selectedInputIndices = {0};
    result.combinedProbability = {.1};
    result.confidence.inputIndices = {1};
    result.confidence.qValues = {.5};
    QVERIFY(CandidateBundleIO::writeReport({loaded}, result, temp.filePath("bad.radiantDIA")) != Error::eNoError);
}

QTEST_MAIN(CandidateBundleIOTests)
#include "CandidateBundleIOTests.moc"
