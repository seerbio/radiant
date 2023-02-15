#include "FastaReader.h"

#include "ErrorUtils.h"
#include "StringUtils.h"

#include<QFile>
#include <QString>
#include <QTextStream>

#include <iostream>


Err FastaReader::parseFastaFile(const QString &filePath, bool addDecoys)
{
    ERR_INIT;

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

            if(!fastaEntry.fastaSequence.isEmpty()){

                fastaEntry.proteinIdNumber = m_fastaEntries.size();
                m_fastaEntries.insert(fastaEntry.proteinIdNumber, fastaEntry);
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
        fastaEntry.proteinIdNumber = m_fastaEntries.size();
        m_fastaEntries.insert(fastaEntry.proteinIdNumber, fastaEntry);
        fastaEntry.clear();
    }

    file.close();

    if (!hasDecoys() && addDecoys) {
        addDecoysToFastaEntrties();
    }

    ERR_RETURN;
}


QMap<ProteinId, FastaEntry> FastaReader::fastaEntries() const {
    return m_fastaEntries;
}


bool FastaReader::hasDecoys() {

    int targetFastaEntryCount = 0;
    int decoyFastaEntrtryCount = 0;

    for (const FastaEntry &fe : m_fastaEntries) {

        if (fe.isDecoy()) {
            decoyFastaEntrtryCount++;
        }
        else {
            targetFastaEntryCount++;
        }
    }

    return decoyFastaEntrtryCount == targetFastaEntryCount;
}


namespace {

    QMap<ProteinId, FastaEntry> filterAnyDecoysThatMightBePresent(const QMap<ProteinId, FastaEntry> &fastaEntries) {

        QMap<ProteinId, FastaEntry> newfastaEntries;

        for (const FastaEntry &fe : fastaEntries) {

            if(fe.isDecoy()) {
                continue;
            }

            newfastaEntries.insert(newfastaEntries.size(), fe);
        }

        return newfastaEntries;
    }


    QString updateFastaDescription(const QString &fastaDescription) {

        const QChar splitter = '|';
        QStringList fastaDescriptionSplit = fastaDescription.split(splitter);
        const QString spOrTr = fastaDescriptionSplit.value(0).replace(">", "");
        const QString newSpOrTr = ">rev_" + spOrTr;
        fastaDescriptionSplit[0] = newSpOrTr;

        return fastaDescriptionSplit.join(splitter);
    }



}//namespace
Err FastaReader::addDecoysToFastaEntrties() {

    ERR_INIT

    const QMap<ProteinId, FastaEntry> newfastaEntries = filterAnyDecoysThatMightBePresent(m_fastaEntries);

    m_fastaEntries = newfastaEntries;

    for (const FastaEntry &fe : newfastaEntries) {

        e = ErrorUtils::isNotEmpty(fe.fastaDescription); ree;
        e = ErrorUtils::isNotEmpty(fe.fastaSequence); ree;

        const QString newFastaDescription = updateFastaDescription(fe.fastaDescription);

        QString ogSequenceToReverse = fe.fastaSequence;
        std::reverse(ogSequenceToReverse.begin(), ogSequenceToReverse.end());

        FastaEntry feRev;
        feRev.fastaDescription = newFastaDescription;
        feRev.fastaSequence = ogSequenceToReverse;
        feRev.proteinIdNumber = m_fastaEntries.size();

        m_fastaEntries.insert( feRev.proteinIdNumber, feRev);
    }

    e = ErrorUtils::isTrue(hasDecoys());

    ERR_RETURN
}
