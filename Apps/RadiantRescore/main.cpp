#include "CandidateBundleIO.h"
#include "CandidateFeatureSchema.h"

#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QSet>

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    if (argc != 2) {
        qCritical() << "Usage: RadiantRescore cohort.json";
        return 2;
    }
    QFile file(QString::fromLocal8Bit(argv[1]));
    if (!file.open(QIODevice::ReadOnly)) return 3;
    const auto configurationBytes = file.readAll();
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(configurationBytes, &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) return 3;
    const auto config = document.object();
    const QSet<QString> allowed{"runs", "output"};
    for (const auto &key : config.keys()) if (!allowed.contains(key)) return 3;
    if (!config["runs"].isArray() || config["runs"].toArray().isEmpty()
        || config["output"].toString().isEmpty() || QFile::exists(config["output"].toString())) return 3;
    QVector<QVector<CandidateBundleIO::View>> views;
    QVector<CandidatePoolRescorer::Run> runs;
    QSet<QString> names;
    QSet<QString> runSourceHashes;
    QString libraryHash, fastaHash;
    for (const auto &value : config["runs"].toArray()) {
        if (!value.isObject()) return 3;
        const auto runConfig = value.toObject();
        for (const auto &key : runConfig.keys()) if (key != "name" && key != "views") return 3;
        const QString name = runConfig["name"].toString();
        if (name.isEmpty() || names.contains(name) || !runConfig["views"].isArray()
            || runConfig["views"].toArray().isEmpty()) return 3;
        names.insert(name);
        CandidatePoolRescorer::Run run;
        run.name = name;
        QVector<CandidateBundleIO::View> runViews;
        QString sourceHash;
        QSet<QString> directories;
        for (const auto &path : runConfig["views"].toArray()) {
            if (!path.isString() || path.toString().isEmpty()) return 3;
            CandidateBundleIO::View view;
            if (CandidateBundleIO::read(path.toString(), &view) != Error::eNoError) {
                qCritical() << "Invalid candidate bundle" << path.toString();
                return 4;
            }
            if (directories.contains(view.directory)) return 3;
            directories.insert(view.directory);
            const QString source = view.provenance["source_sha256"].toString();
            const QString library = view.provenance["library_sha256"].toString();
            const QString fasta = view.provenance["fasta_sha256"].toString();
            if (source.size() != 64 || library.size() != 64 || fasta.size() != 64) return 4;
            if (sourceHash.isEmpty()) sourceHash = source;
            if (libraryHash.isEmpty()) libraryHash = library;
            if (fastaHash.isEmpty()) fastaHash = fasta;
            if (source != sourceHash || library != libraryHash || fasta != fastaHash) {
                qCritical() << "Views must share one raw source per run and one library/FASTA across runs";
                return 4;
            }
            if (qint64(run.identities.size()) + view.candidates.identities.size() > INT_MAX) return 4;
            const int viewNumber = runViews.size();
            for (auto id : view.candidates.identities) {
                id.view = viewNumber;
                run.identities.push_back(std::move(id));
            }
            run.evidence.append(view.candidates.evidence);
            run.rawFeatures.append(view.candidates.rawFeatures);
            runViews.push_back(std::move(view));
        }
        if (runSourceHashes.contains(sourceHash)) {
            qCritical() << "The same raw source cannot appear as multiple runs";
            return 4;
        }
        runSourceHashes.insert(sourceHash);
        runs.push_back(std::move(run));
        views.push_back(std::move(runViews));
    }
    CandidatePoolRescorer::Settings settings;
    QVector<CandidatePoolRescorer::Result> results;
    if (CandidatePoolRescorer::score(runs, settings, &results) != Error::eNoError) {
        qCritical() << "Combined candidate scoring failed";
        return 5;
    }
    const QDir output(config["output"].toString());
    const QFileInfo outputInfo(output.absolutePath());
    QDir parent(outputInfo.absolutePath());
    if (!QDir().mkpath(parent.absolutePath()) || !parent.mkdir(outputInfo.fileName())) return 6;
    QJsonArray summaries;
    for (int index = 0; index < results.size(); ++index) {
        const auto &result = results[index];
        QJsonObject summary{{"name", result.name}, {"selected_candidates", result.selectedInputIndices.size()},
                            {"retained_precursors_including_decoys", result.confidence.inputIndices.size()}};
        QJsonArray inputViews;
        for (const auto &view : views[index]) {
            inputViews.push_back(QJsonObject{
                {"directory", view.directory}, {"manifest_sha256", view.manifestSha256},
                {"files_sha256", view.fileSha256}, {"rows", view.candidates.identities.size()},
                {"source_sha256", view.provenance["source_sha256"]}});
        }
        summary["input_views"] = inputViews;
        int targets = 0, decoys = 0;
        for (int row = 0; row < result.confidence.inputIndices.size(); ++row) {
            if (result.confidence.qValues[row] <= .01) {
                runs[index].identities[result.confidence.inputIndices[row]].isDecoy ? ++decoys : ++targets;
            }
        }
        summary["targets_at_precursor_q_0_01"] = targets;
        summary["decoys_at_precursor_q_0_01"] = decoys;
        if (!result.confidence.inputIndices.isEmpty()) {
            const QString name = QString::number(index) + ".radiantDIA";
            if (CandidateBundleIO::writeReport(
                views[index], result, output.filePath(name), settings.neuralNet.folds) != Error::eNoError)
                return 6;
            summary["report"] = name;
            summary["report_sha256"] = CandidateBundleIO::fileHash(output.filePath(name));
        } else {
            summary["report"] = QJsonValue::Null;
        }
        summaries.push_back(summary);
    }
    const auto &nn = settings.neuralNet;
    const QJsonObject neuralNet{
        {"folds", nn.folds}, {"networks", nn.networks}, {"epochs", nn.epochs},
        {"threads", nn.threads}, {"seed", nn.seed}, {"learning_rate", nn.learningRate},
        {"nodes_fraction", nn.nodesFraction}, {"focal_loss_gamma", nn.focalLossGamma},
        {"log_intensities", nn.logIntensities}, {"shuffle_each_epoch", nn.shuffleEachEpoch},
        {"non_tims_features", PeptideFamilyNeuralNet::nonTimsFeatures().size()},
        {"normalization", "Min-max across each selected model pool"},
        {"family_hash", "Unmodified origin; I to L; SHA256 first8 hex"}};
    QJsonObject software{
        {"executable_sha256", CandidateBundleIO::fileHash(QCoreApplication::applicationFilePath())}};
    const QString runtimeManifest = CandidateBundleIO::fileHash(
        QDir(QCoreApplication::applicationDirPath()).filePath("manifest.json"));
    if (!runtimeManifest.isEmpty()) software["runtime_manifest_sha256"] = runtimeManifest;
    const QJsonObject manifest{
        {"format", "radiant-combined-candidate-results"}, {"version", 1},
        {"library_sha256", libraryHash}, {"fasta_sha256", fastaHash},
        {"configuration", config}, {"runs", summaries},
        {"configuration_sha256", QString::fromLatin1(
            QCryptographicHash::hash(configurationBytes, QCryptographicHash::Sha256).toHex())},
        {"software", software}, {"feature_schema_sha256", CandidateFeatureSchema::id()},
        {"neural_net", neuralNet}, {"minimum_shared_fragments", settings.minimumSharedFragments},
        {"original_candidate_limit", settings.originalCandidateLimit},
        {"expanded_candidate_limit", settings.expandedCandidateLimit},
        {"combination", "Equal-logit original pooled/per-run, then equal-logit with expanded pooled on original prefix"},
        {"confidence", "Per-run precursor only; fresh exact-tie +1 q after peak and physical selection"},
        {"peptide_protein_confidence", "Unassigned (q=1, best flags=0); require a separate downstream analysis"},
        {"decoy_ratio_diagnostic", "Unassigned (-1); exported LDA confidence is not reused"},
        {"match_between_runs", false}};
    QFile completion(output.filePath("manifest.json"));
    if (!completion.open(QIODevice::WriteOnly | QIODevice::NewOnly)) return 6;
    const auto bytes = QJsonDocument(manifest).toJson(QJsonDocument::Indented);
    if (completion.write(bytes) != bytes.size() || !completion.flush()) return 6;
    qInfo() << "Wrote combined precursor results to" << output.absolutePath();
    return 0;
}
