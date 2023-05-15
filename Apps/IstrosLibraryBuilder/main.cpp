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

    const QString model1FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_1.hdf5.json");
    const QString model2FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_2.hdf5.json");
    const QString model3FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_3.hdf5.json");
    const QString model4FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_4.hdf5.json");

    QString fragLibFilePath;

    LibraryBuilderWorkFlow libraryBuilderWorkFlow;
    e = libraryBuilderWorkFlow.init(
            pythiaParameters,
            model1FilePath,
            model2FilePath,
            model3FilePath,
            model4FilePath
    );
    if (e != eNoError) {
        qDebug() << "buiding Library builder input succeeded";
    }

//    e = libraryBuilderWorkFlow.exec(
//            outputPepLibBuilderCSVFilePath,
//            &fragLibFilePath
//    );
//    if (e != eNoError) {
//        qDebug() << "Building library failed";
//    }

    return 0;

}


