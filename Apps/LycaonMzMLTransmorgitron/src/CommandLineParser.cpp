//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_MZML_PATH = QStringLiteral("mzML-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_MZML_PATH, QObject::tr("*.mzML file"));
}

bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 2) {
        qCritical() << QStringLiteral("Expected 1 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.mzMLFilePath = args[1];

    const bool mzMLPathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.mzMLFilePath,
            {S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION}
    );

    if (!mzMLPathIsValid) {
        qCritical() << QStringLiteral("First command line argument *.mzML, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_MZML_PATH  << m_cliParams.mzMLFilePath;

    return m_cliParams;
}
