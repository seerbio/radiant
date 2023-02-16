//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "Error.h"
#include "workflows/PythiaFraggerWorkflow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QDirIterator>
#include <QElapsedTimer>


using namespace Error;


Err buildPythiaParameters(const QString &pythiaParametersFilePath,
                          PythiaParameters *pythiaParameters) {

    ERR_INIT

    PythiaParameterReader pythiaParameterReader;
    e = pythiaParameterReader.readFile(pythiaParametersFilePath); ree;

    e = pythiaParameterReader.loadPythiaParameters(pythiaParameters); ree;
    pythiaParameters->print();

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

    PythiaParameters pythiaParameters;
    e = buildPythiaParameters(cliParameters.pythiaParametersFilePath, &pythiaParameters);
    if (e != eNoError) {
        qDebug() << QStringLiteral("Unable to load parameters from %1")
        .arg(cliParameters.pythiaParametersFilePath);
        return 1;
    }

    PythiaFraggerWorkflow pythiaFraggerWorkflow;
    e = pythiaFraggerWorkflow.init(pythiaParameters,
                                   cliParameters.fastaFilePath);
    if (e != eNoError) {
        qCritical() << e;
        qDebug() << "Unable to init workflow";
        return 1;
    }

    QFileInfo fileInfo(cliParameters.dataFilePath);
    const bool isFile = fileInfo.isFile();
    const bool isDirectory = fileInfo.isDir();

    if (isFile) {

        e = pythiaFraggerWorkflow.processDataFile(cliParameters.dataFilePath, true);
        if (e != eNoError) {
            qCritical() << e;
            return 1;
        }

        qDebug() << "All files processed in" << et.elapsed() / 1000 << "seconds";
        return 0;
    }

    else if (isDirectory) {

        QDirIterator it(cliParameters.dataFilePath);
        while (it.hasNext()) {

            const QString &dataFilePath = it.next();
            const QFileInfo fi(dataFilePath);
            const QString fileExtension = fi.suffix();

            if (StringUtils::stringsMatch(FileReadersLibNameSpace::MZML_FILE_EXTENSION, fileExtension, false) ||
                StringUtils::stringsMatch(FileReadersLibNameSpace::HDF_FILE_EXTENSION, fileExtension, false)) {

                e = pythiaFraggerWorkflow.processDataFile(dataFilePath, true);
                if (e != eNoError) {
                    qDebug() << "*** IONS2 FILE NOT FRAGGED PROPERLY ****";
                    qDebug() << "ions2FilePath:" << dataFilePath;
                    qDebug() << "*****************************************";
                    eee_absorb;
                }
            }
        }

        qDebug() << "All files processed in" << et.elapsed() / 1000 << "seconds";
        return 0;
    }

    return 1;
}


