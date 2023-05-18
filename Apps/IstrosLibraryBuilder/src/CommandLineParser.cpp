//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_PEPTIDES_CSV_PATH = QStringLiteral("peptides-csv-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_PEPTIDES_CSV_PATH, QObject::tr("*.csv file"));
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

    m_cliParams.peptideLibBuilderCSVFilePath = args[1];

    const bool fastaFilePathIsValid = CommandLineParserUtils::checkFileNameExtension(
            m_cliParams.peptideLibBuilderCSVFilePath,
            S_GLOBAL_SETTINGS.CSV_FILE_EXTENSION
            );

    if (!fastaFilePathIsValid) {
        qCritical() << QStringLiteral("First command line argument *.fasta, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[2];

    const bool pythiaPathIsValid = CommandLineParserUtils::checkFileNameExtension(
            m_cliParams.pythiaParametersFilePath,
            S_GLOBAL_SETTINGS.PYTHIA_FILE_EXTENSION
            );

    if (!pythiaPathIsValid) {
        qCritical() << QStringLiteral("Second command line argument *.pythia argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_PEPTIDES_CSV_PATH << m_cliParams.peptideLibBuilderCSVFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;

    return m_cliParams;
}
