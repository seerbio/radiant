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

    const QString fragLibPath = cliParameters.fragLibFilePath;
    const QString fastaFilePath = cliParameters.fastaFilePath;
    const QString pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QString msDataFile = cliParameters.msDataFile;

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
            fastaFilePath
    );
    if (e != eNoError) {
        qDebug() << "Error initializing Pythia Workflow Libraries";
        return 1;
    }

    // e = pythiaDiaFFWorkflow.processFile(msDataFile);
    // if (e != eNoError) {
    //     qDebug() << msDataFile << "Did not run completely";
    //     return 1;
    // }

    qDebug() << "PSMing done in" << et.elapsed() << "mSec";

    return 0;
}
