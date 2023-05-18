//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "PythiaDIAWorkflow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QDirIterator>
#include <QElapsedTimer>


using namespace Error;


Err buildPythiaParameters(
        const QString &pythiaParametersFilePath,
        PythiaParameters *pythiaParameters
        ) {

    ERR_INIT

    PythiaParameterReader pythiaParameterReader;
    e = pythiaParameterReader.readFile(pythiaParametersFilePath); ree;

    e = pythiaParameterReader.loadPythiaParameters(pythiaParameters); ree;
    pythiaParameters->print();

    ERR_RETURN
}

Err getDataFilesFromDirectory(
        const QString &dataFilesDirectory,
        QStringList *dataFiles
        ) {

    ERR_INIT

    QDirIterator it(dataFilesDirectory);
    while (it.hasNext()) {

        const QString &dataFilePath = it.next();
        const QFileInfo fi(dataFilePath);
        const QString fileExtension = fi.suffix();

        const QString prqSuffix = S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION;

        if (StringUtils::stringsMatch(prqSuffix, fileExtension, false)) {
            dataFiles->push_back(dataFilePath);
        }
    }

    ERR_RETURN
}

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
    const QString pythiaParamsFilePath = cliParameters.pythiaParametersFilePath;
    const QString msDataFilesDirectory = cliParameters.msDataFilesDirectory;

    PythiaParameters pythiaParameters;
    e = buildPythiaParameters(
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
            fragLibPath
    );
    if (e != eNoError) {
        qDebug() << "Error initializing Pythia Workflow Libraries";
        return 1;
    }

    QStringList dataFiles;
    e = getDataFilesFromDirectory(
            msDataFilesDirectory,
            &dataFiles
            ); ree;

    e = ErrorUtils::isNotEmpty(dataFiles); ree;

    for (const QString &dataFilePath : dataFiles) {

        e = pythiaDiaWorkflow.processFile(dataFilePath);
        if (e != eNoError) {
            qDebug() << dataFilePath << "Did not run completely";
            return 1;
        }
    }

    qDebug() << "PSMing done in" << et.elapsed() << "mSec";

    return 0;

}
