//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"

#include "ParallelUtils.h"

#include <QElapsedTimer>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeCoor, int> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


FragLibReader::FragLibReader()
: m_mzMin(200.0)
, m_mzMax(2000.0) {}

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
        e = ErrorUtils::isEqual(row.mzVals.size(), row.intensityVals.size()); ree;

        const QStringList ionLabels = row.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

        e = ErrorUtils::isEqual(row.mzVals.size(), ionLabels.size()); ree;

        for (int i = 0; i < row.mzVals.size(); i++) {
            MS2Ion ion;
            ion.mz = row.mzVals.at(i);
            ion.intensity = row.intensityVals.at(i);
            ion.ionLabel = ionLabels.at(i);
            ion.iRT = static_cast<float>(row.iRT);
            ms2Ions->push_back(ion);
        }

        ERR_RETURN
    }

}//namespace
Err FragLibReader::getMS2Ions(
        double massStart,
        double massEnd,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
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

    const int maxMs2Ions = 12; //TODO change this
    e = fragLibReaderRowsToMs2IonsMap(
            fragLibReaderRows,
            maxMs2Ions,
            peptideSequenceChargeKeyVsCandidatePeptide
            ); ree;

    qDebug() << "MS2 Predictions count:" << peptideSequenceChargeKeyVsCandidatePeptide->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}

Err FragLibReader::getMS2Ions(QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide) {

    ERR_INIT

    const int maxMs2Ions = 300;
    e = getMS2IonsTopN(
            maxMs2Ions,
            peptideSequenceChargeKeyVsCandidatePeptide
            ); ree;

    ERR_RETURN
}

Err FragLibReader::getMS2IonsTopN(
        int topNMs2Ions,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = ParquetReader::read(
            m_fragLibFilePath,
            &fragLibReaderRows
    ); ree

    e = fragLibReaderRowsToMs2IonsMap(
            fragLibReaderRows,
            topNMs2Ions,
            peptideSequenceChargeKeyVsCandidatePeptide
    ); ree

    qDebug() << "All MS2 Predictions count:" << peptideSequenceChargeKeyVsCandidatePeptide->size() << "retrieved in" << et.elapsed() << "mSec";

    ERR_RETURN
}

void FragLibReader::filterMs2IonsByMz(
        double mzStart,
        double mzEnd,
        QVector<MS2Ion> *ms2Ions
) {

    const auto terminatorLogic = [mzStart, mzEnd](const MS2Ion &ion){
        return !(mzStart <= ion.mz && ion.intensity <= mzEnd);
    };

    const auto terminator = std::remove_if(
            ms2Ions->begin(),
            ms2Ions->end(),
            terminatorLogic
    );

    ms2Ions->erase(terminator, ms2Ions->end());
}

namespace {

    void removeMzVals(
            double mzMin,
            double mzMax,
            QVector<MS2Ion> *ms2Ions
            ) {

        const auto terminatorLogic = [mzMin, mzMax](const MS2Ion &ms2Ion){
            return ms2Ion.mz < mzMin || ms2Ion.mz > mzMax;
        };

        const auto terminator = std::remove_if(
                ms2Ions->begin(),
                ms2Ions->end(),
                terminatorLogic
                );

        ms2Ions->erase(terminator, ms2Ions->end());
    }

}//namespace
void FragLibReader::getTopNMostIntenseMs2Ions(
        int topNMs2Ions,
        double mzMin,
        double mzMax,
        QVector<MS2Ion> *ms2Ions
) {

    removeMzVals(mzMin, mzMax, ms2Ions);

    const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
        return l.intensity < r.intensity;
    };

    std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);

    topNMs2Ions = std::min(topNMs2Ions, ms2Ions->size());

    ms2Ions->resize(topNMs2Ions);

    for (int i = 0; i < ms2Ions->size(); i++) {
        (*ms2Ions)[i].rank = i;
    }

}

void FragLibReader::filterMs2IonsByIntensity(
        double intensityThreshold,
        QVector<MS2Ion> *ms2Ions
        ) {

    const auto terminatorLogic = [intensityThreshold](const MS2Ion &ion){
        return ion.intensity < intensityThreshold;
    };

    const auto terminator = std::remove_if(
            ms2Ions->begin(),
            ms2Ions->end(),
            terminatorLogic
            );

    ms2Ions->erase(terminator, ms2Ions->end());

}

namespace {

    double calculateMassFromThomson(
            double mz,
            int charge
    ) {
        const double PROTON = 1.0072;
        return (mz * charge) - (charge * PROTON);
    }

}//namespace
Err FragLibReader::buildFragIonLibForCandidates(
        const QString &fragLibUri,
        int chargeMin,
        int chargeMax,
        double mzTargetMin,
        double mzTargetMax,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
        ) {

    ERR_INIT

    FragLibReader fragLibReader;
    e = fragLibReader.init(fragLibUri); ree;

    for (Charge charge = chargeMin; charge <= chargeMax; ++charge) {

        const double massStart = calculateMassFromThomson(
                mzTargetMin,
                charge
        );

        const double massEnd = calculateMassFromThomson(
                mzTargetMax,
                charge
        );

        e = fragLibReader.getMS2Ions(
                massStart,
                massEnd,
                peptideSequenceChargeKeyVsCandidatePeptide
        ); ree

    }

    ERR_RETURN
}

Err FragLibReader::peptideStringWithModsFromPeptideSequenceChargeKey(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        PeptideStringWithMods *peptideStringWithMods,
        Charge *charge
){

    ERR_INIT

    const int expectedSplitSize = 2;

    const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
            S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
            QString::SkipEmptyParts
    );

    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize); ree;

    *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
    ); ree

    ERR_RETURN
}

Err FragLibReader::fragLibReaderRowsToMs2IonsMap(
        const QVector<FragLibReaderRow> &fragLibReaderRows,
        int topNMs2Ions,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
) const {

    ERR_INIT

    for (const FragLibReaderRow &row : fragLibReaderRows) {

        e = ErrorUtils::isFalse(peptideSequenceChargeKeyVsCandidatePeptide->contains(row.peptideSequenceChargeKey)); ree;

        QVector<MS2Ion> ms2Ions;
        e = buildMzIons(
                row,
                &ms2Ions
        ); ree;

        CandidatePeptide candidatePeptide;
        candidatePeptide.totalFragmentCount = ms2Ions.size();

        getTopNMostIntenseMs2Ions(
                topNMs2Ions,
                m_mzMin,
                m_mzMax,
                &ms2Ions
                ); ree

        candidatePeptide.ms2Ions = ms2Ions;

        candidatePeptide.mass = row.mass;
        candidatePeptide.iRt = row.iRT;
        candidatePeptide.isDecoy = row.isDecoy;

        peptideSequenceChargeKeyVsCandidatePeptide->insert(row.peptideSequenceChargeKey, candidatePeptide);

    }

    ERR_RETURN
}

namespace {

    RTree loadRtree(const QMap<PeptideString, CandidatePeptide> &peptideSequenceVsMS2Ions) {

        std::vector<rTreePoint> cloudLoader;

        for (auto it = peptideSequenceVsMS2Ions.begin(); it != peptideSequenceVsMS2Ions.end(); it++) {

            const PeptideStringWithMods &peptideStr = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value().ms2Ions;

            for (const MS2Ion &ms2Ion : ms2Ions) {
                rTreeCoor coor(ms2Ion.mz, 0.0);
                cloudLoader.emplace_back(coor, -1);
            }
        }

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    QMap<MzHashed, MZION> buildMzHashedVsMzIon(const QList<CandidatePeptide> &candidatePeptides) {

        QMap<MzHashed, MZION> mzHashedVsMzIon;

        const int precision = S_GLOBAL_SETTINGS.HASHING_PRECISION;

        for (const CandidatePeptide &candidatePeptide : candidatePeptides) {

            for (const MS2Ion &ms2Ion : candidatePeptide.ms2Ions) {

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, precision);
                const MZION mzIon = MathUtils::pRound(ms2Ion.mz, precision);

                mzHashedVsMzIon.insert(mzHashed, mzIon);
            }
        }

        return mzHashedVsMzIon;
    }

}
Err FragLibReader::generateFragmentFrequencies(
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        double ppmTol,
        QMap<MzHashed, FrequencyPercent> *fragmentFrequencies
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(peptideStringWithModsVsCandidatePeptide); ree
    fragmentFrequencies->clear();

    RTree rTree = loadRtree(peptideStringWithModsVsCandidatePeptide);
    const QMap<MzHashed, MZION> mzIonHashedVsMzIon
        = buildMzHashedVsMzIon(peptideStringWithModsVsCandidatePeptide.values());

    const int rTreeIonCount = rTree.size();
    e = ErrorUtils::isNotEqual(rTreeIonCount, 0); ree

    for (auto it = mzIonHashedVsMzIon.begin(); it != mzIonHashedVsMzIon.end(); it++) {

        const MzHashed mzHashed = it.key();
        const MZION mzIon = it.value();

        const double mzTol = MathUtils::calculatePPM(mzIon, ppmTol);
        const double mzMin = mzIon - mzTol;
        const double mzMax = mzIon + mzTol;

        const rTreeSearchBox queryBox(
                rTreeCoor(mzMin, 0.0),
                rTreeCoor(mzMax, 0.0)
        );

        std::vector<rTreePoint> result;
        rTree.query(bgi::intersects(queryBox), std::back_inserter(result));

        const double fragmentFreq = static_cast<double>(result.size()) / rTreeIonCount;

        fragmentFrequencies->insert(mzHashed, fragmentFreq);
    }

    qDebug() << "Frequences calculated in" << et.elapsed() << "mSec";

    ERR_RETURN
}
