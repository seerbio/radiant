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


struct FILEREADERSLIB_EXPORTS FastaEntry{

   QString fastaDescription;
   QString fastaSequence;
   QString proteinId;
   QString geneName;

   int proteinIdNumber = -1;

   [[nodiscard]] bool isDecoy() const {
       return fastaDescription.contains("rev_", Qt::CaseInsensitive)
            || fastaDescription.contains("decoy_", Qt::CaseInsensitive);
   }

   void clear(){
       fastaDescription.clear();
       fastaSequence.clear();
   }

    friend QDataStream &operator <<(QDataStream &stream, const FastaEntry &fe) {
        stream << fe.fastaDescription;
        stream << fe.fastaSequence;
        stream << fe.proteinId;
        stream << fe.geneName;
        stream << fe.proteinIdNumber;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, FastaEntry &fe) {
        stream >> fe.fastaDescription;
        stream >> fe.fastaSequence;
        stream >> fe.proteinId;
        stream >> fe.geneName;
        stream >> fe.proteinIdNumber;

        return stream;
    }

};


class FILEREADERSLIB_EXPORTS FastaReader
{
public:

    FastaReader() = default;
    ~FastaReader() = default;

    Err parseFastaFile(const QString &filePath, bool addDecoys);

    bool hasDecoys();

    [[nodiscard]] QMap<ProteinId, FastaEntry> fastaEntries() const;

private:

    Err addDecoysToFastaEntrties();

private:
    QMap<ProteinId, FastaEntry> m_fastaEntries;

};

#endif // FASTAREADER_H
