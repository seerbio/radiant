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
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("config-path");
    const QString ARG_DATAFILE_PATH = QStringLiteral("data-file-path");
    const QString ARG_OUTPUT_FOLDER = QStringLiteral("output-folder");
} // END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_FRAGLIB_PATH, QObject::tr("*.fragLibFF file"));
    addPositionalArgument(ARG_FASTA_PATH, QObject::tr("*.fasta file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.radiantConfig file"));
    addPositionalArgument(ARG_DATAFILE_PATH, QObject::tr("data directory path"));

    addOption(QCommandLineOption(
        {ARG_OUTPUT_FOLDER},
        QObject::tr("Set the output folder path."),
        QObject::tr("path")
    ));
}

bool CommandLineParser::validateArguments(const QStringList &args) {
    // Process the arguments first to avoid "call process() or parse() before isSet" error
    this->process(args);

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
            {
                S_GLOBAL_SETTINGS.DOT_FASTA.mid(1, -1),
                S_GLOBAL_SETTINGS.DOT_FAS.mid(1, -1)
            }
    );
    if (!fragLibBackgroundFilePathIsValid) {
        qCritical() << QStringLiteral("Second command line argument: expected a FASTA file with .fasta or .fas extension");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[3];
    const bool pythiaPathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.pythiaParametersFilePath,
            {
                S_GLOBAL_SETTINGS.CONFIG_FILE_EXTENSION,
                S_GLOBAL_SETTINGS.LEGACY_CONFIG_FILE_EXTENSION,
                S_GLOBAL_SETTINGS.TOML_FILE_EXTENSION,
            }
            );
    if (!pythiaPathIsValid) {
        qCritical() << QStringLiteral("Third command line argument *.radiantConfig argument invalid");
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

    if (isSet(ARG_OUTPUT_FOLDER)) {
        m_cliParams.outputFolderPath = value(ARG_OUTPUT_FOLDER);
        m_cliParams.outputFolderPath = m_cliParams.outputFolderPath.back() == "/"
                                     ? m_cliParams.outputFolderPath
                                     : m_cliParams.outputFolderPath + "/";

        QDir dir;
        if (!dir.exists(m_cliParams.outputFolderPath)) {
            if (dir.mkpath(m_cliParams.outputFolderPath)) {
                qDebug() << "Directory created successfully:" << m_cliParams.outputFolderPath;
            } else {
                qDebug() << "Failed to create directory:" << m_cliParams.outputFolderPath;
            }
        }
    }

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_FRAGLIB_PATH << m_cliParams.fragLibFilePath;
    qDebug() << ARG_FASTA_PATH << m_cliParams.fastaFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;
    qDebug() << ARG_DATAFILE_PATH  << m_cliParams.msDataFile;
    qDebug() << ARG_OUTPUT_FOLDER  << m_cliParams.outputFolderPath;

    return m_cliParams;
}