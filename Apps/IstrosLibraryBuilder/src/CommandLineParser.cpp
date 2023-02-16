//
// Created by Drucifer on 11/24/2021.
//

#include "CommandLineParser.h"

#include "FileReadersLibNameSpace.h"
#include "StringUtils.h"

#include <QDebug>
#include <QFileInfo>

using namespace FileReadersLibNameSpace;

namespace {
    const QString ARG_RAW_DATA_PATH = QStringLiteral("raw-data-path");
    const QString ARG_PYTHIA_PARAMS = QStringLiteral("pythia-path");
    const QString ARG_FASTA = QStringLiteral("fasta-path");
}//END NAMESPACE


CommandLineParser::CommandLineParser() {
    addHelpOption();
    addPositionalArgument(ARG_RAW_DATA_PATH, QObject::tr("*.mzml file"));
    addPositionalArgument(ARG_PYTHIA_PARAMS, QObject::tr("*.pythia file"));
    addPositionalArgument(ARG_FASTA, QObject::tr("*.fasta file"));
}


namespace {

    bool checkFileNameExtension(const QString &filePath,
                                const QString &expectedFileExtension) {

        QFileInfo fi(filePath);
        const QString fileSuffix = fi.suffix();

        if (StringUtils::stringsMatch(fileSuffix, expectedFileExtension, false) && fi.isFile()) {
            return true;
        }

        else if (fi.isDir()) {
            return true;
        }

        return false;
    }


}//namespace
bool CommandLineParser::validateArguments(const QStringList &args) {

    QStringList argumentsLocal(args);
    if (argumentsLocal.size() != 4) {
        qCritical() << QStringLiteral("Expected 3 arguments.  Received %1").arg(argumentsLocal.size() - 1);
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.dataFilePath = args[1];
    const bool mzmlPathOrDirIsValid
        = checkFileNameExtension(m_cliParams.dataFilePath, MZML_FILE_EXTENSION);

    const bool hdfPathOrDirIsValid
            = checkFileNameExtension(m_cliParams.dataFilePath, HDF_FILE_EXTENSION);


    if (!(mzmlPathOrDirIsValid || hdfPathOrDirIsValid)) {
        qCritical() << QStringLiteral("First command line argument *.mzml, *.hdf / directory argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.pythiaParametersFilePath = args[2];
    const bool pythiaPathIsValid
            = checkFileNameExtension(m_cliParams.pythiaParametersFilePath, PYTHIA_FILE_EXTENSION);
    if (!pythiaPathIsValid) {
        qCritical() << QStringLiteral("Second command line argument *.pythia argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    m_cliParams.fastaFilePath = args[3];
    const bool fastaPathIsValid
            = checkFileNameExtension(m_cliParams.fastaFilePath, FASTA_FILE_EXTENSION);
    if (!fastaPathIsValid) {
        qCritical() << QStringLiteral("Thrid command line argument *.fasta argument invalid");
        argumentsLocal.append("-h");
        process(argumentsLocal);
        return false;
    }

    process(argumentsLocal);

    return true;
}


const CommandLineParser::CliParameters &CommandLineParser::getCliParams() const {

    qDebug() << "IONS2 PATH" << m_cliParams.dataFilePath;
    qDebug() << "PYTHIA PARAMS PATH"  << m_cliParams.pythiaParametersFilePath;
    qDebug() << "FASTA PATH" << m_cliParams.fastaFilePath;

    return m_cliParams;
}
