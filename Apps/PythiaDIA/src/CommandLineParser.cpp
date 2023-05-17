//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_PEPTIDE_CSV_DATA_PATH = QStringLiteral("peptides-csv-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_PEPTIDE_CSV_DATA_PATH, QObject::tr("*.csv file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.pythia file"));
}


bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 3) {
        qCritical() << QStringLiteral("Expected 3 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.peptidesCSVFilePath = args[1];

//    const bool peptidesCSVFilePathIsValid
//        = checkFileNameExtension(m_cliParams.peptidesCSVFilePath, S_GLOBAL_SETTINGS.DOT_CSV);
//
//    if (!peptidesCSVFilePathIsValid) {
//        qCritical() << QStringLiteral("First command line argument *.csv, argument invalid");
//        argumentsLocal.append("-h");
//        process(argumentsLocal);
//        return false;
//    }

    m_cliParams.pythiaParametersFilePath = args[2];

//    const bool pythiaPathIsValid
//            = checkFileNameExtension(m_cliParams.pythiaParametersFilePath, ".pythia");
//    if (!pythiaPathIsValid) {
//        qCritical() << QStringLiteral("Second command line argument *.pythia argument invalid");
//        argumentsLocal.append("-h");
//        process(argumentsLocal);
//        return false;
//    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << "PEPTIDES PATH" << m_cliParams.peptidesCSVFilePath;
    qDebug() << "PYTHIA PARAMS PATH"  << m_cliParams.pythiaParametersFilePath;

    return m_cliParams;
}
