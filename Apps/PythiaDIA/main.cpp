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

    const QString mzMLFileURI
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    const QString fragLibPath
            = "/home/anichols/Desktop/2022_02_22_Homo_sapiens_UP000005640.fragLib";

    PythiaDIAWorkflow pythiaDiaWorkflow;
    e = pythiaDiaWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            fragLibPath
    );
    if (e != eNoError) {
        qDebug() << "you done messed up";
        //TODO properly handle error
    }

    e = pythiaDiaWorkflow.processFile(mzMLFileURI);
    if (e != eNoError) {
        qDebug() << "you done messed up";
        //TODO properly handle error
    }


    qDebug() << "PSMing done in" << et.elapsed() << "mSec";

    return 0;

}


