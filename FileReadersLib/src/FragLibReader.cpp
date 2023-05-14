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
            QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions,
            QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy
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

            peptideSequenceChargeKeyVsIsDecoy->insert(
                    row.peptideSequenceChargeKey,
                    row.isDecoy
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
        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideSequenceChargeKeyVsMS2Ions,
        QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy
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
            peptideSequenceChargeKeyVsMS2Ions,
            peptideSequenceChargeKeyVsIsDecoy
            ); ree;

    qDebug() << peptideSequenceChargeKeyVsMS2Ions->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}


void FragLibReader::filterMs2IonsByMz(
        double mzStart,
        double mzEnd,
        QVector<MS2Ion> *ms2Ions
) {

    const auto terminatorLogic = [mzStart, mzEnd](const MS2Ion &ion){
        return !(mzStart <= ion.x() && ion.x() <= mzEnd);
    };

    const auto terminator = std::remove_if(
            ms2Ions->begin(),
            ms2Ions->end(),
            terminatorLogic
    );

    ms2Ions->erase(terminator, ms2Ions->end());
}

void FragLibReader::getTopNMostIntenseMs2Ions(
        int topNMs2Ions,
        QVector<MS2Ion> *ms2Ions
) {

    const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
        return l.y() < r.y();
    };

    std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);

    topNMs2Ions = std::min(topNMs2Ions, ms2Ions->size());

    ms2Ions->resize(topNMs2Ions);

    const auto sortMzAsc = [](const MS2Ion &l, const MS2Ion &r) {
        return l.x() < r.x();
    };

    std::sort(ms2Ions->begin(), ms2Ions->end(), sortMzAsc);
}