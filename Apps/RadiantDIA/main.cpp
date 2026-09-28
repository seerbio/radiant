//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "src/CandidateCohort.h"
#include "CommandLineParserUtils.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "PythiaDIAFFWorkflow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QtConcurrent>
#include <QThreadPool>

using namespace Error;


int main(int argc, char *argv[]) {

    ERR_INIT

//    QThreadPool::globalInstance()->setMaxThreadCount(8); //TODO make this settable.

    QElapsedTimer et;
    et.start();

    QCoreApplication app(argc, argv);
    const auto arguments = QCoreApplication::arguments();
    if (arguments.size() > 1 && arguments[1] == "--candidate-cohort") {
        if (arguments.size() != 3) {
            qCritical() << "Usage: RadiantDIA --candidate-cohort cohort.json";
            return 2;
        }
        return runCandidateCohort(arguments[2]);
    }
    CommandLineParser parser;

    if (!parser.validateArguments(QCoreApplication::arguments())) {
        return 1;
    }

    const CommandLineParser::CliParameters &cliParameters = parser.getCliParams();

    const QString &fragLibPath = cliParameters.fragLibFilePath;
    const QString &fastaFilePath = cliParameters.fastaFilePath;
    const QString &pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QString &msDataFile = cliParameters.msDataFile;
    const QString &outputFolderPath = cliParameters.outputFolderPath;

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

    e = pythiaDiaFFWorkflow.processFile(msDataFile);
    if (e != eNoError) {
        qDebug() << msDataFile << "Did not run completely";
        return 1;
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMing done in" << et.elapsed() << "mSec";

    return 0;
}
