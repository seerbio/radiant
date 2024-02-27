//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"
#include "FragLibTsvReader.h"

#include <QElapsedTimer>
#include <QFileInfo>

Err FragLibReader::getFragLibReaderRows(
        const QString &fragLibFilePath,
        double massStart,
        double massEnd,
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(massEnd > massStart); ree;

    QFileInfo fi(fragLibFilePath);
    const QString fileSuffix = fi.suffix();

    const QStringList viableFileExtensions = {
            FragLibReaderNamespace::FRAG_LIB_FF_SUFFIX.toLower(),
            FragLibReaderNamespace::FRAG_LIB_CSV_FF_SUFFIX.toLower(),
            FragLibReaderNamespace::FRAG_LIB_TSV_FF_SUFFIX.toLower()
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
                FragLibReaderNamespace::MASS,
                {massStart, massEnd},
                fragLibReaderRows
        ); ree;

        ERR_RETURN
    }

    FragLibTsvReader fragLibTsvReader;
    e = fragLibTsvReader.getFragLibReaderRows(
            fragLibFilePath,
            massStart,
            massEnd,
            fragLibReaderRows
            ); ree;

    qDebug() << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}
