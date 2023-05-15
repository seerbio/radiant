//
// Created by Drucifer on 12/31/2021.
//

#include "src/CommandLineParser.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "Error.h"
#include "LibraryBuilderWorkFlow.h"
#include "PythiaParameterReader.h"
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

Err generatePeptidesFromFasta(
        const PythiaParameters &pythiaParameters,
        const QString &fastaFilePath,
        const QString &targetMzCollisionCSVFilePath,
        const QString &outputFilePath
        ) {

    ERR_INIT

    FastaFileToPeptidesListWorkFlow fastaFileToPeptidesListWorkFlow;

    e = fastaFileToPeptidesListWorkFlow.init(pythiaParameters); ree;

    e = fastaFileToPeptidesListWorkFlow.exec(
            fastaFilePath,
            targetMzCollisionCSVFilePath,
            outputFilePath
            ); ree;

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
    e = buildPythiaParameters(
            cliParameters.pythiaParametersFilePath,
            &pythiaParameters
            );
    if (e != eNoError) {
        qDebug() << "Error reading Pythia Parameters file!";
    }

    const QString outputPepLibBuilderCSVFilePath
        = cliParameters.fastaFilePath + S_GLOBAL_SETTINGS.DOT_CSV;

    e = generatePeptidesFromFasta(
            pythiaParameters,
            cliParameters.fastaFilePath,
            cliParameters.targetMzCollisionCSVFilePath,
            outputPepLibBuilderCSVFilePath
            );
    if (e != eNoError) {
        qDebug() << "Error generating peptides from fasta";
    }

    return 0;

}
