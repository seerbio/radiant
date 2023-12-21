//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"

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

    e = ErrorUtils::isTrue(
            FragLibReaderNamespace::FRAG_LIB_FF_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
            ); ree;

    fragLibReaderRows->clear();

    QElapsedTimer et;
    et.start();

    e = ParquetReader::read(
            fragLibFilePath,
            FragLibReaderNamespace::MASS,
            {massStart, massEnd},
            fragLibReaderRows
            ); ree;

    qDebug() << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}
