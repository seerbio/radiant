//
// Created by Drucifer on 11/24/2021.
//
#include "CommandLineParser.h"

#include "CommandLineParserUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <QFileInfo>

namespace {
    const QString ARG_LIBRARY_PATH = QStringLiteral("library-path");
    const QString ARG_FASTA_PATH = QStringLiteral("fasta-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("config-path");
    const QString ARG_DATAFILE_PATH = QStringLiteral("data-file-path");
    const QString ARG_OUTPUT_FOLDER = QStringLiteral("output-folder");
    const QString ARG_THERMO_RAW_READER = QStringLiteral("thermo-raw-reader");
} // END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_LIBRARY_PATH, QObject::tr("library file (e.g. *.tsv)"));
    addPositionalArgument(ARG_FASTA_PATH, QObject::tr("*.fasta file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.radiantConfig file"));
    addPositionalArgument(ARG_DATAFILE_PATH, QObject::tr("data path (*.mzML, Thermo *.raw, or Bruker *.d)"));

    addOption(QCommandLineOption(
        {ARG_OUTPUT_FOLDER},
        QObject::tr("Set the output folder path."),
        QObject::tr("path")
    ));
    addOption(QCommandLineOption(
        {ARG_THERMO_RAW_READER},
        QObject::tr("RadiantThermoRawReader executable or assembly path for direct Thermo RAW input "
                    "(default: RADIANT_THERMO_RAW_READER or PATH)."),
        QObject::tr("path")
    ));
}

bool CommandLineParser::validateArguments(const QStringList &args) {
    // Process the arguments first to avoid "call process() or parse() before isSet" error
    this->process(args);

    QStringList argumentsLocal(args);
    const QStringList positional = positionalArguments();
    if (positional.size() < 4) {
        qCritical() << QStringLiteral("Expected 4 arguments.  Received %1").arg(positional.size());
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fragLibFilePath = positional[0];
    const bool fragLibFilePathIsValid = CommandLineParserUtils::checkFileNameExtensions(
            m_cliParams.fragLibFilePath,
            {
                S_GLOBAL_SETTINGS.DOT_FRAGLIB_FF.mid(1, 9),
                S_GLOBAL_SETTINGS.DOT_TSV.mid(1,3),
                S_GLOBAL_SETTINGS.DOT_SPECLIB.mid(1,7)
            }
            );
    if (!fragLibFilePathIsValid) {
        qCritical() << QStringLiteral("First command line argument must be a library file (e.g. *.tsv)");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fastaFilePath = positional[1];
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

    m_cliParams.pythiaParametersFilePath = positional[2];
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

    m_cliParams.msDataFile = positional[3];
    const bool dataFilesPathIsValid = CommandLineParserUtils::isMassSpectrometryDataPath(m_cliParams.msDataFile);
    if (!dataFilesPathIsValid) {
        qCritical() << QStringLiteral("Fourth command line argument must be *.mzML, Thermo *.raw, or Bruker *.d");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.thermoRawReaderPath = value(ARG_THERMO_RAW_READER);

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

    qDebug() << ARG_LIBRARY_PATH << m_cliParams.fragLibFilePath;
    qDebug() << ARG_FASTA_PATH << m_cliParams.fastaFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;
    qDebug() << ARG_DATAFILE_PATH  << m_cliParams.msDataFile;
    qDebug() << ARG_OUTPUT_FOLDER  << m_cliParams.outputFolderPath;
    if (!m_cliParams.thermoRawReaderPath.isEmpty()) {
        qDebug() << ARG_THERMO_RAW_READER << m_cliParams.thermoRawReaderPath;
    }

    return m_cliParams;
}
