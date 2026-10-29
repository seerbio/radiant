//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "CommandLineParserUtils.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "PythiaDIAFFWorkflow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QSet>
#include <QtConcurrent>
#include <QThreadPool>

using namespace Error;


int main(int argc, char *argv[]) {

    ERR_INIT

//    QThreadPool::globalInstance()->setMaxThreadCount(8); //TODO make this settable.

    QElapsedTimer et;
    et.start();

    QCoreApplication app(argc, argv);
    CommandLineParser parser;

    if (!parser.validateArguments(QCoreApplication::arguments())) {
        return 1;
    }

    const CommandLineParser::CliParameters &cliParameters = parser.getCliParams();

    if (!cliParameters.thermoRawReaderPath.isEmpty()) {
        qputenv("RADIANT_THERMO_RAW_READER", cliParameters.thermoRawReaderPath.toUtf8());
    }

    const QString &fragLibPath = cliParameters.fragLibFilePath;
    const QString &fastaFilePath = cliParameters.fastaFilePath;
    const QString &pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QStringList &msDataFiles = cliParameters.msDataFiles;
    const QString &outputFolderPath = cliParameters.outputFolderPath;

    // Reject colliding destinations before processing so a batch cannot silently overwrite itself.
    QSet<QString> outputPaths;
    for (const QString &msDataFile : msDataFiles) {
        const QFileInfo inputInfo(QDir::cleanPath(msDataFile));
        if (!CommandLineParserUtils::isMassSpectrometryDataPath(msDataFile)
            || !inputInfo.exists()) {
            continue;
        }
        const QDir destination = outputFolderPath.isEmpty()
            ? QDir(inputInfo.dir().canonicalPath())
            : QDir(QDir(outputFolderPath).canonicalPath());
        const QString outputPath = destination.filePath(
            inputInfo.fileName() + S_GLOBAL_SETTINGS.DOT_RADIANT_DIA_FILE_EXTENSION);
        if (outputPaths.contains(outputPath)) {
            qCritical().noquote() << "Batch inputs would overwrite the same output:" << outputPath;
            return 1;
        }
        outputPaths.insert(outputPath);
    }

    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
            pythiaParamsFilePath,
            &pythiaParameters
            );
    if (e != eNoError) {
        qDebug() << "Error reading pythia parameters";
        return 1;
    }

    PythiaDIAFFWorkflow pythiaDiaFFWorkflow;
    e = pythiaDiaFFWorkflow.init(
            pythiaParameters,
            fragLibPath,
            fastaFilePath,
            outputFolderPath
    );
    if (e != eNoError) {
        qDebug() << "Error initializing Pythia Workflow Libraries";
        return 1;
    }

    int completedCount = 0;
    int failedCount = 0;
    for (int index = 0; index < msDataFiles.size(); ++index) {
        const QString &msDataFile = msDataFiles[index];
        qInfo().noquote() << QStringLiteral("Processing file %1/%2: %3")
                            .arg(index + 1).arg(msDataFiles.size()).arg(msDataFile);
        if (!CommandLineParserUtils::isMassSpectrometryDataPath(msDataFile)) {
            qCritical().noquote() << "FAILED:" << msDataFile
                                 << "- expected mzML, Thermo raw, Bruker d, or prqFF input";
            ++failedCount;
            continue;
        }
        e = pythiaDiaFFWorkflow.processFile(msDataFile);
        if (e != eNoError) {
            qCritical().noquote() << "FAILED:" << msDataFile << "- error" << static_cast<int>(e);
            ++failedCount;
            continue;
        }
        qInfo().noquote() << "COMPLETED:" << msDataFile;
        ++completedCount;
    }

    qInfo().noquote() << QStringLiteral("Batch summary: %1 completed, %2 failed, %3 total")
                        .arg(completedCount).arg(failedCount).arg(msDataFiles.size());
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMing done in" << et.elapsed() << "mSec";

    return failedCount == 0 ? 0 : 1;
}
