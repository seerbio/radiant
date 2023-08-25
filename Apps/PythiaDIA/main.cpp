//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "CommandLineParserUtils.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "PythiaDIAWorkflow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QElapsedTimer>


using namespace Error;


int main(int argc, char *argv[]) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    QCoreApplication app(argc, argv);
    CommandLineParser parser;

    if (!parser.validateArguments(QCoreApplication::arguments())) {
        return 1;
    }

    const CommandLineParser::CliParameters &cliParameters = parser.getCliParams();

    const QString fragLibPath = cliParameters.fragLibFilePath;
    const QString fragLibBackgroundPath = cliParameters.fragLibFilePath;
    const QString pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QString msDataFilesDirectory = cliParameters.msDataFilesDirectory;

    //TODO make sure this optional argument works.
    const QString iRTReCalFilePath = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/TestData.iRT");


    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
            pythiaParamsFilePath,
            &pythiaParameters
            );
    if (e != eNoError) {
        qDebug() << "Error reading pythia parameters";
        return 1;
    }

    PythiaDIAWorkflow pythiaDiaWorkflow;
    e = pythiaDiaWorkflow.init(
            pythiaParameters,
            fragLibPath,
            iRTReCalFilePath
    );
    if (e != eNoError) {
        qDebug() << "Error initializing Pythia Workflow Libraries";
        return 1;
    }

//    QStringList dataFiles;
//    e = CommandLineParserUtils::getDataFilesFromDirectory(
//            msDataFilesDirectory,
//            &dataFiles
//            );
//    if (e != eNoError) {
//        qDebug() << "Error reading data files.";
//        return 1;
//    }

//    e = ErrorUtils::isNotEmpty(dataFiles);
//    if (e != eNoError) {
//        qDebug() << "No data files found.";
//        return 1;
//    }

    qDebug() << "DKFJDSL" << cliParameters.msDataFilesDirectory;
    e = pythiaDiaWorkflow.processFile(cliParameters.msDataFilesDirectory);
//    if (e != eNoError) {
//        qDebug() << cliParameters.msDataFilesDirectory << "Did not run completely";
//        return 1;
//    }

    qDebug() << "PSMing done in" << et.elapsed() << "mSec";

    return 0;

}
