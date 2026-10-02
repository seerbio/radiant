//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_MSDATA_PATH = QStringLiteral("ms-data-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_MSDATA_PATH, QObject::tr("*.mzML file or Bruker *.d directory"));
}

bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 2) {
        qCritical() << QStringLiteral("Expected 1 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.msDataFilePath = args[1];

    const bool msDataPathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.msDataFilePath,
            {
                S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION,
                S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION
            }
    );

    if (!msDataPathIsValid) {
        qCritical() << QStringLiteral("First command line argument must be *.mzML or a Bruker *.d directory");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_MSDATA_PATH  << m_cliParams.msDataFilePath;

    return m_cliParams;
}
