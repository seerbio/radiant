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

    /**
    * @brief Parses a FASTA file in FastaReader.
    *
    * This function parses a FASTA file specified by the `filePath` parameter in the FastaReader
    * instance. It checks if the file exists, opens it for reading, and then reads and processes
    * each line of the file to extract FASTA entries. The extracted entries are stored in the
    * `m_fastaEntries` member variable. The function ensures that unwanted characters are removed
    * from the sequence and that entries are only added if they are not decoy sequences. Any
    * encountered errors during the process are indicated by the returned Err code.
    *
    * @param filePath The file path of the FASTA file to be parsed.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err parseFastaFile(const QString &filePath);

    [[nodiscard]] QMap<ProteinId, FastaEntry> fastaEntries() const;


private:
    QMap<ProteinId, FastaEntry> m_fastaEntries;

};

#endif // FASTAREADER_H
