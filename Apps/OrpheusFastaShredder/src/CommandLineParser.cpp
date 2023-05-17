//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_FASTA_PATH = QStringLiteral("fasta-path");
    const QString ARG_TARGET_CE_CSV_PATH = QStringLiteral("target-mz-ce-csv-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_FASTA_PATH, QObject::tr("*.fasta file"));
    addPositionalArgument(ARG_TARGET_CE_CSV_PATH, QObject::tr("*.csv file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.pythia file"));
}

bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 4) {
        qCritical() << QStringLiteral("Expected 3 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fastaFilePath = args[1];

    const bool fastaFilePathIsValid
        = CommandLineParserUtils::checkFileNameExtension(m_cliParams.fastaFilePath, "fasta");

    if (!fastaFilePathIsValid) {
        qCritical() << QStringLiteral("First command line argument *.fasta, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.targetMzCollisionCSVFilePath = args[2];

    const bool targetMzCollisionCSVFilePathValid
            = CommandLineParserUtils::checkFileNameExtension(m_cliParams.targetMzCollisionCSVFilePath, "csv");
    if (!targetMzCollisionCSVFilePathValid) {
        qCritical() << QStringLiteral("Second command line argument *.csv argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[3];

    const bool pythiaPathIsValid
            = CommandLineParserUtils::checkFileNameExtension(m_cliParams.pythiaParametersFilePath, "pythia");
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

    qDebug() << ARG_FASTA_PATH << m_cliParams.fastaFilePath;
    qDebug() << ARG_TARGET_CE_CSV_PATH << m_cliParams.targetMzCollisionCSVFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;

    return m_cliParams;
}
