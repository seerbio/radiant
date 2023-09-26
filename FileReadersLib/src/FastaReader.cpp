#include "FastaReader.h"

#include "ErrorUtils.h"
#include "StringUtils.h"

#include<QFile>
#include <QString>
#include <QTextStream>


namespace {
    bool isDecoy(const QString &fastaDescription) {
        return fastaDescription.contains("rev_", Qt::CaseInsensitive)
               || fastaDescription.contains("decoy_", Qt::CaseInsensitive);
    }
}//namespace
Err FastaReader::parseFastaFile(const QString &filePath)
{
    ERR_INIT;

    e = ErrorUtils::fileExists(filePath); ree;

    QFile file(filePath);
    e = ErrorUtils::isTrue(file.open(QIODevice::ReadOnly), eFileError); ree;

    const QString ALLOWED_FASTA_CHARS = QStringLiteral("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

    QTextStream in(&file);

    FastaEntry fastaEntry;

    while(!in.atEnd()) {

        QString line = in.readLine().trimmed();
        if(line.isEmpty()){
            continue;
        }

        if(line.at(0) == '>'){

            if(!fastaEntry.fastaSequence.isEmpty() && !isDecoy(fastaEntry.fastaSequence)
            ){
                m_fastaEntries.insert(m_fastaEntries.size(), fastaEntry);
                fastaEntry.clear();
            }

            fastaEntry.fastaDescription = line;
        }

        else{
            line = line.toUpper();
            StringUtils::removeUnwantedChars(ALLOWED_FASTA_CHARS, &line);
            fastaEntry.fastaSequence += line;
        }
    }

    if(in.atEnd()){
        m_fastaEntries.insert(m_fastaEntries.size(), fastaEntry);
        fastaEntry.clear();
    }

    file.close();

    ERR_RETURN;
}


QMap<ProteinId, FastaEntry> FastaReader::fastaEntries() const {
    return m_fastaEntries;
}
