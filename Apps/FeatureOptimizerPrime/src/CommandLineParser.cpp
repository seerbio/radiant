//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>


namespace {
    const QString ARG_FRAGLIB_PATH = QStringLiteral("fraglibff-path");
    const QString ARG_FASTA_PATH = QStringLiteral("fasta-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
    const QString ARG_DATAFILE_PATH = QStringLiteral("data-file-path");
}//END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_FRAGLIB_PATH, QObject::tr("*.fragLibFF file"));
    addPositionalArgument(ARG_FASTA_PATH, QObject::tr("*.fasta file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.pythia file"));
    addPositionalArgument(ARG_DATAFILE_PATH, QObject::tr("data directory path"));
}


bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() < 5) {
        qCritical() << QStringLiteral("Expected 4 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fragLibFilePath = args[1];
    const bool fragLibFilePathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.fragLibFilePath,
            {
                S_GLOBAL_SETTINGS.DOT_FRAGLIB_FF.mid(1, 9),
                S_GLOBAL_SETTINGS.DOT_CSV.mid(1,3),
                S_GLOBAL_SETTINGS.DOT_TSV.mid(1,3),
                S_GLOBAL_SETTINGS.DOT_SPECLIB.mid(1,7)
            }
            );
    if (!fragLibFilePathIsValid) {
        qCritical() << QStringLiteral("First command line argument *.fragLibFF, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fastaFilePath = args[2];
    const bool fragLibBackgroundFilePathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.fastaFilePath,
            {S_GLOBAL_SETTINGS.DOT_FASTA.mid(1, 7)}
    );
    if (!fragLibBackgroundFilePathIsValid) {
        qCritical() << QStringLiteral("Second command line argument *.fragLib, argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[3];
    const bool pythiaPathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.pythiaParametersFilePath,
            {S_GLOBAL_SETTINGS.CONFIG_FILE_EXTENSION}
            );
    if (!pythiaPathIsValid) {
        qCritical() << QStringLiteral("Third command line argument *.pythiaConfig argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.msDataFile = args[4];
    const bool dataFilesPathIsValid = CommandLineParserUtils::checkFileNameExtensions(m_cliParams.msDataFile, {"prqFF"})
                                   || CommandLineParserUtils::checkFileNameExtensions(m_cliParams.msDataFile, {"mzML"})
                                   || CommandLineParserUtils::checkFileNameExtensions(m_cliParams.msDataFile, {"d"});
    if (!dataFilesPathIsValid) {
        qCritical() << QStringLiteral("Fourth command line argument data directory path argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_FRAGLIB_PATH << m_cliParams.fragLibFilePath;
    qDebug() << ARG_FASTA_PATH << m_cliParams.fastaFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;
    qDebug() << ARG_DATAFILE_PATH  << m_cliParams.msDataFile;

    return m_cliParams;
}
