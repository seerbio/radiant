//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "CommandLineParserUtils.h"
#include "Error.h"
#include "ConvertMzMLToParquetWorkFlow.h"
#include "StringUtils.h"

#include <QCoreApplication>
#include <QDirIterator>
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

    QString outputFilePath;
    e = ConvertMzMLToParquetWorkFlow::convertMzMLToParquet(
            cliParameters.msDataFilePath,
            &outputFilePath
            );
    if (e != eNoError) {
        qDebug() << outputFilePath << "did not convert correctly";
        return 1;
    }

    qDebug() << outputFilePath << "convered in" << et.elapsed() << "mSec";

    return 0;
}
