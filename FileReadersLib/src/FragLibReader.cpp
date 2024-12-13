//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"
#include "FragLibTsvReader.h"
#include "SpecLibReader.h"

#include <QElapsedTimer>
#include <QFileInfo>

Err FragLibReader::getFragLibReaderRows(
        const QString &fragLibFilePath,
        QVector<FragLibReaderRow> *fragLibReaderRows
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

    if (fragLibFilePath.contains(S_GLOBAL_SETTINGS.DOT_FRAGLIB_FF)) {
        e = ParquetReader::read(
                fragLibFilePath,
                fragLibReaderRows
        ); ree;
    }

    else if (fragLibFilePath.contains(S_GLOBAL_SETTINGS.DOT_SPECLIB)) {
        e = SpecLibReader::getFragLibReaerRows(
                fragLibFilePath,
                fragLibReaderRows
        ); ree;
    }

    // FragLibTsvReader fragLibTsvReader;
    // e = fragLibTsvReader.getFragLibReaderRows(
    //         fragLibFilePath,
    //         fragLibReaderRows
    //         ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}
