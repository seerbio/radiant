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
    const QString ARG_PER_FILE_LOGS = QStringLiteral("per-file-logs");
} // END NAMESPACE

CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_LIBRARY_PATH, QObject::tr("library file (e.g. *.tsv)"));
    addPositionalArgument(ARG_FASTA_PATH, QObject::tr("*.fasta file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.radiantConfig file"));
    addPositionalArgument(ARG_DATAFILE_PATH,
                          QObject::tr("One or more data paths (*.mzML, Thermo *.raw, or Bruker *.d), "
                                      "processed sequentially with one library load."),
                          QStringLiteral("<data-file-path> [data-file-path ...]"));

    addOption(QCommandLineOption(
        {ARG_OUTPUT_FOLDER},
        QObject::tr("Set the output folder path."),
        QObject::tr("path")
    ));
    addOption(QCommandLineOption(
        {ARG_PER_FILE_LOGS},
        QObject::tr("Also write each input's processing diagnostics to "
                    "<input-filename>.log alongside its output.")
    ));
    QCommandLineOption thermoRawReaderOption(
        {ARG_THERMO_RAW_READER},
        QObject::tr("RadiantThermoRawReader executable or assembly path for direct Thermo RAW input "
                    "(default: RADIANT_THERMO_RAW_READER or PATH)."),
        QObject::tr("path")
    );
    thermoRawReaderOption.setFlags(QCommandLineOption::HiddenFromHelp);
    addOption(thermoRawReaderOption);
}

bool CommandLineParser::validateArguments(const QStringList &args) {
    // Process the arguments first to avoid "call process() or parse() before isSet" error
    this->process(args);

    const QStringList positional = positionalArguments();
    if (positional.size() < 4) {
        qCritical() << QStringLiteral("Expected library, FASTA, config, and at least one data path. "
                                     "Received %1 arguments").arg(positional.size());
        qInfo().noquote() << helpText();
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
        qInfo().noquote() << helpText();
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
        qInfo().noquote() << helpText();
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
        qInfo().noquote() << helpText();
        return false;
    }

    m_cliParams.msDataFiles = positional.mid(3);

    m_cliParams.thermoRawReaderPath = value(ARG_THERMO_RAW_READER);
    m_cliParams.perFileLogs = isSet(ARG_PER_FILE_LOGS);

    if (isSet(ARG_OUTPUT_FOLDER)) {
        m_cliParams.outputFolderPath = value(ARG_OUTPUT_FOLDER);
        if (m_cliParams.outputFolderPath.isEmpty()) {
            qCritical() << "The output folder path must not be empty";
            return false;
        }
        m_cliParams.outputFolderPath = QDir::cleanPath(m_cliParams.outputFolderPath) + "/";

        QDir dir;
        if (!dir.exists(m_cliParams.outputFolderPath)) {
            if (dir.mkpath(m_cliParams.outputFolderPath)) {
                qDebug() << "Directory created successfully:" << m_cliParams.outputFolderPath;
            } else {
                qCritical() << "Failed to create directory:" << m_cliParams.outputFolderPath;
                return false;
            }
        }
    }

    return true;
}

const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << ARG_LIBRARY_PATH << m_cliParams.fragLibFilePath;
    qDebug() << ARG_FASTA_PATH << m_cliParams.fastaFilePath;
    qDebug() << ARG_PYTHIA_PARAMS  << m_cliParams.pythiaParametersFilePath;
    qDebug() << ARG_DATAFILE_PATH  << m_cliParams.msDataFiles;
    qDebug() << ARG_OUTPUT_FOLDER  << m_cliParams.outputFolderPath;
    if (!m_cliParams.thermoRawReaderPath.isEmpty()) {
        qDebug() << ARG_THERMO_RAW_READER << m_cliParams.thermoRawReaderPath;
    }

    return m_cliParams;
}
