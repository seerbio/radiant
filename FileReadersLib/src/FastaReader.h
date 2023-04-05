#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"

#include <QDataStream>
#include <QMap>
#include <QString>
#include <QVector>


using namespace Error;


struct FILEREADERSLIB_EXPORTS FastaEntry {

   QString fastaDescription;
   QString fastaSequence;
   QString proteinId;
   QString geneName;

   void clear(){
       fastaDescription.clear();
       fastaSequence.clear();
   }

};


class FILEREADERSLIB_EXPORTS FastaReader
{
public:

    FastaReader() = default;
    ~FastaReader() = default;

    Err parseFastaFile(const QString &filePath);

    [[nodiscard]] QMap<ProteinId, FastaEntry> fastaEntries() const;


private:
    QMap<ProteinId, FastaEntry> m_fastaEntries;

};

#endif // FASTAREADER_H
