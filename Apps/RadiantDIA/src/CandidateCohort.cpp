#include "CandidateCohort.h"

#include "PythiaDIAFFWorkflow.h"
#include "PythiaParameterReader.h"

#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSet>

namespace {
bool keysAllowed(const QJsonObject &object, const QSet<QString> &allowed) {
    for (const auto &key : object.keys())
        if (!allowed.contains(key)) return false;
    return true;
}
}

int runCandidateCohort(const QString &configurationPath) {
    QElapsedTimer timer;
    timer.start();
    QFile file(configurationPath);
    if (!file.open(QIODevice::ReadOnly)) return 3;
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(file.readAll(), &error);
    if (error.error != QJsonParseError::NoError || !document.isObject()) return 3;
    const auto configuration = document.object();
    if (!keysAllowed(configuration, {"library", "fasta", "config", "runs"})) return 3;
    const QDir base = QFileInfo(configurationPath).absoluteDir();
    const auto path = [&](const QJsonValue &value) -> QString {
        if (!value.isString() || value.toString().isEmpty()) return {};
        return QDir::cleanPath(base.absoluteFilePath(value.toString()));
    };
    const QString libraryPath = path(configuration["library"]);
    const QString fastaPath = path(configuration["fasta"]);
    const QString parameterPath = path(configuration["config"]);
    if (!QFileInfo(libraryPath).isFile() || !QFileInfo(fastaPath).isFile()
        || !QFileInfo(parameterPath).isFile()) return 3;
    PythiaParameters parameters;
    if (PythiaParameterReader::buildPythiaParameters(parameterPath, &parameters) != eNoError
        || !parameters.isValid() || !parameters.candidateBundleOnly
        || parameters.candidateBundleLimit < 2) return 3;
    if (!configuration["runs"].isArray() || configuration["runs"].toArray().isEmpty()) return 3;
    struct Run {
        QString input;
        QVector<PythiaDIAFFWorkflow::CandidateView> views;
    };
    QVector<Run> runs;
    QSet<QString> inputs, outputs;
    for (const auto &entry : configuration["runs"].toArray()) {
        if (!entry.isObject()) return 3;
        const auto object = entry.toObject();
        if (!keysAllowed(object, {"input", "views"})) return 3;
        Run run;
        run.input = path(object["input"]);
        const QFileInfo input(run.input);
        if (!input.isFile() || inputs.contains(input.canonicalFilePath())) return 3;
        inputs.insert(input.canonicalFilePath());
        if (!object["views"].isArray() || object["views"].toArray().isEmpty()) return 3;
        for (const auto &value : object["views"].toArray()) {
            if (!value.isObject()) return 3;
            const auto view = value.toObject();
            if (!keysAllowed(view, {"minimum_fragments", "shared_fragments", "output"})) return 3;
            PythiaDIAFFWorkflow::CandidateView parsed;
            parsed.minimumFragments = view["minimum_fragments"].toInt(-1);
            parsed.sharedFragments = view["shared_fragments"].toInt(-1);
            parsed.outputDirectory = path(view["output"]);
            if (!view["minimum_fragments"].isDouble() || !view["shared_fragments"].isDouble()
                || view["minimum_fragments"].toDouble() != parsed.minimumFragments
                || view["shared_fragments"].toDouble() != parsed.sharedFragments
                || parsed.minimumFragments < 3 || parsed.minimumFragments > 12
                || parsed.sharedFragments < 2 || parsed.sharedFragments > 12
                || parsed.outputDirectory.isEmpty() || outputs.contains(parsed.outputDirectory)
                || QFileInfo::exists(parsed.outputDirectory)) return 3;
            outputs.insert(parsed.outputDirectory);
            run.views.push_back(std::move(parsed));
        }
        runs.push_back(std::move(run));
    }
    PythiaDIAFFWorkflow::LibraryHandle library;
    if (PythiaDIAFFWorkflow::prepareLibrary(
            libraryPath, parameters.useAlternativeDecoys, &library) != eNoError) return 4;
    for (const auto &run : runs) {
        // Each raw file gets fresh mutable search state. Only parsed library
        // rows are shared across files; calibrated data stay local to one file.
        PythiaDIAFFWorkflow workflow;
        if (workflow.init(parameters, libraryPath, fastaPath, {}, library) != eNoError) return 4;
        if (workflow.processCandidateViews(run.input, run.views) != eNoError) return 5;
    }
    qInfo() << "Candidate cohort complete:" << runs.size() << "files in" << timer.elapsed() << "msec";
    return 0;
}
