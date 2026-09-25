//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"
#include "AminoAcids.h"
#include "FragLibTsvReader.h"
#include "PeptideStringWithMods.h"
#include "SpecLibReader.h"

#include <QElapsedTimer>
#include <QFileInfo>
#include <algorithm>

namespace {

    bool hasValidSequence(const FragLibReaderRow &fragLibReaderRow) {
        const QStringList keySplit = fragLibReaderRow.peptideSequenceChargeKey.split(
            S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
            Qt::SkipEmptyParts
            );

        if (keySplit.size() != 2) {
            return false;
        }

        const PeptideString sequenceNoMods = PeptideStringWithMods(keySplit.front()).removeUniModChars();
        return AminoAcids::validPeptideSequence(sequenceNoMods);
    }

    void filterInvalidSequenceEntries(QList<FragLibReaderRow> *fragLibReaderRows) {
        const auto terminator = std::remove_if(
            fragLibReaderRows->begin(),
            fragLibReaderRows->end(),
            [](const FragLibReaderRow &fragLibReaderRow) {
                return !hasValidSequence(fragLibReaderRow);
            }
            );
        fragLibReaderRows->erase(terminator, fragLibReaderRows->end());
    }

}

Err FragLibReader::getFragLibReaderRows(
        const QString &fragLibFilePath,
        bool useAlternativeDecoys,
        QList<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;

    QFileInfo fi(fragLibFilePath);
    const QString fileSuffix = fi.suffix();

    const QStringList viableFileExtensions = {
            FragLibReaderNamespace::FRAG_LIB_FF_SUFFIX.toLower(),
            FragLibReaderNamespace::FRAG_LIB_CSV_FF_SUFFIX.toLower(),
            FragLibReaderNamespace::FRAG_LIB_TSV_FF_SUFFIX.toLower(),
            FragLibReaderNamespace::SPEC_LIB_CSV_FF_SUFFIX.toLower()
    };

    e = ErrorUtils::isTrue(
            viableFileExtensions.contains(fileSuffix.toLower()),
            eFileIncorrectTypeError
            ); ree;

    fragLibReaderRows->clear();

    QElapsedTimer et;
    et.start();

    if (fragLibFilePath.endsWith(S_GLOBAL_SETTINGS.DOT_FRAGLIB_FF)) {
        e = ParquetReader::read(
                fragLibFilePath,
                fragLibReaderRows
        ); ree;
    }

    else if (fragLibFilePath.endsWith(S_GLOBAL_SETTINGS.DOT_SPECLIB)) {
        e = SpecLibReader::getFragLibReaerRows(
                fragLibFilePath,
                fragLibReaderRows
        ); ree;
    }
    else if (fragLibFilePath.endsWith(S_GLOBAL_SETTINGS.DOT_TSV)){
            FragLibTsvReader fragLibTsvReader;
            fragLibTsvReader.setEnableTerminalByPenultimateDecoyAnnotationShift(useAlternativeDecoys);
            e = fragLibTsvReader.getFragLibReaderRows(
                    fragLibFilePath,
                    fragLibReaderRows
                    ); ree;
    }

    const int countBeforeFilter = fragLibReaderRows->size();
    filterInvalidSequenceEntries(fragLibReaderRows);
    const int invalidSequenceCount = countBeforeFilter - fragLibReaderRows->size();
    if (invalidSequenceCount > 0) {
        qWarning() << qPrintable(S_GLOBAL_TIMER.elapsed())
                   << "Dropped"
                   << invalidSequenceCount
                   << "library entries with invalid peptide sequences";
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}
