//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"

#include "ParallelUtils.h"

#include <QElapsedTimer>

Err FragLibReader::init(const QString &fragLibFilePath) {
    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    m_fragLibFilePath = fragLibFilePath;

    ERR_RETURN
}

namespace {

    Err buildMzIons(
            const FragLibReaderRow &row,
            QVector<MS2Ion> *ms2Ions
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(row.mzVals); ree;
        e = ErrorUtils::isEqual(
                row.mzVals.size(),
                row.intensityVals.size()
                ); ree;

        for (int i = 0; i < row.mzVals.size(); i++) {
            ms2Ions->push_back({row.mzVals.at(i), row.intensityVals.at(i)});
        }

        ERR_RETURN
    }

    Err fragLibReaderRowsToMs2IonsMap(
            const QVector<FragLibReaderRow> &fragLibReaderRows,
            QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions
            ) {

        ERR_INIT

        for (const FragLibReaderRow &row : fragLibReaderRows) {

            QVector<MS2Ion> ms2Ions;
            e = buildMzIons(
                    row,
                    &ms2Ions
                    ); ree;

            e = ErrorUtils::isFalse(peptideStringWithModsVsMS2Ions->contains(row.peptideSequenceChargeKey)); ree;

            peptideStringWithModsVsMS2Ions->insert(
                    row.peptideSequenceChargeKey,
                    ms2Ions
                    );

        }

        ERR_RETURN
    }

//    void filterMs2IonsByNMostIntense(QVector<FragLibReaderRow> *fragLibReaderRows) {
//
//        const auto terminatorLogic = [](const FragLibReader &row){
//            return
//        };
//
//    }

}//namespace
Err FragLibReader::getMS2Ions(
        double massStart,
        double massEnd,
        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideSequenceChargeKeyVsMS2Ions
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = ParquetReader::read(
            m_fragLibFilePath,
            FragLibReaderNamespace::MASS,
            {massStart, massEnd},
            &fragLibReaderRows
            ); ree;

    e = fragLibReaderRowsToMs2IonsMap(
            fragLibReaderRows,
            peptideSequenceChargeKeyVsMS2Ions
            ); ree;

    qDebug() << peptideSequenceChargeKeyVsMS2Ions->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}
