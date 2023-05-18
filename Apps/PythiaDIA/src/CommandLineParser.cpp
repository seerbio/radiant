//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_FRAGLIB_PATH = QStringLiteral("fraglib-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
    const QString ARG_DATAFILES_DIR_PARAMS = QStringLiteral("data-directory-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_FRAGLIB_PATH, QObject::tr("*.fragLib file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.pythia file"));
    addPositionalArgument(ARG_DATAFILES_DIR_PARAMS, QObject::tr("data directory path"));
}


bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 4) {
        qCritical() << QStringLiteral("Expected 3 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fragLibFilePath = args[1];

    const bool fragLibFilePathIsValid = CommandLineParserUtils::checkFileNameExtension(
            m_cliParams.fragLibFilePath,
            S_GLOBAL_SETTINGS.DOT_FRAGLIB.mid(1,7)
            );

    if (!fragLibFilePathIsValid) {
        qCritical() << QStringLiteral("First command line argument *.fragLib, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[2];

    const bool pythiaPathIsValid = CommandLineParserUtils::checkFileNameExtension(
            m_cliParams.pythiaParametersFilePath, S_GLOBAL_SETTINGS.PYTHIA_FILE_EXTENSION);

    if (!pythiaPathIsValid) {
        qCritical() << QStringLiteral("Second command line argument *.pythia argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.msDataFilesDirectory = args[3];

    const bool dataFilesPathIsValid = CommandLineParserUtils::checkFileNameExtension(
            m_cliParams.msDataFilesDirectory,
            ""
            );

    if (!dataFilesPathIsValid) {
        qCritical() << QStringLiteral("Third command line argument data directory path argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_FRAGLIB_PATH << m_cliParams.fragLibFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;
    qDebug() << ARG_DATAFILES_DIR_PARAMS  << m_cliParams.msDataFilesDirectory;

    return m_cliParams;
}
