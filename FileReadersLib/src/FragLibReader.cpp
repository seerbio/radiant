//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"


//#include "MolecularFormula.h"
//#include "ParallelUtils.h"

#include <QElapsedTimer>

//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point.hpp>
//#include <boost/geometry/geometries/box.hpp>
//#include <boost/geometry/index/rtree.hpp>
//
//#include <iostream>
//
//namespace bg = boost::geometry;
//namespace bgi = boost::geometry::index;
//
//using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
//using rTreeSearchBox = bg::model::box<rTreeCoor>;
//using rTreePoint = std::pair<rTreeCoor, int> ;
//using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;



Err FragLibReader::getFragLibReaderRows(
        const QString &fragLibFilePath,
        double massStart,
        double massEnd,
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(massEnd > massStart); ree;

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


//Err FragLibReader::getMS2Ions(QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide) {
//
//    ERR_INIT
//
//    const int maxMs2Ions = 300;
//    e = getMS2IonsTopN(
//            maxMs2Ions,
//            m_mzMin,
//            m_mzMax,
//            m_useBYOnly,
//            peptideSequenceChargeKeyVsCandidatePeptide
//            ); ree;
//
//    ERR_RETURN
//}

//Err FragLibReader::getMS2IonsTopN(
//        int topNMs2Ions,
//        double mzMin,
//        double mzMax,
//        bool byIonsOnly,
//        QVector<FragLibReaderRow> *fragLibReaderRows
//        ) {
//
//    ERR_INIT
//
//    QElapsedTimer et;
//    et.start();
//
//    e = fragLibReaderRowsToMs2IonsMap(
//            m_fragLibReaderRows,
//            m_aminoAcids,
//            topNMs2Ions,
//            mzMin,
//            mzMax,
//            byIonsOnly,
//            peptideSequenceChargeKeyVsCandidatePeptide
//    ); ree
//
//    qDebug() << "All MS2 Predictions count:" << peptideSequenceChargeKeyVsCandidatePeptide->size() << "retrieved in" << et.elapsed() << "mSec";
//
//    ERR_RETURN
//}


//namespace {
//
//    void removeMzVals(
//            double mzMin,
//            double mzMax,
//            QVector<MS2Ion> *ms2Ions
//            ) {
//
//        const auto terminatorLogic = [mzMin, mzMax](const MS2Ion &ms2Ion){
//            return ms2Ion.mz < mzMin || ms2Ion.mz > mzMax;
//        };
//
//        const auto terminator = std::remove_if(
//                ms2Ions->begin(),
//                ms2Ions->end(),
//                terminatorLogic
//                );
//
//        ms2Ions->erase(terminator, ms2Ions->end());
//    }
//
//    void removeNeutralLosses(QVector<MS2Ion> *ms2Ions) {
//
//        const auto terminatorLogic = [](const MS2Ion &ms2Ion){
//            return ms2Ion.ionLabel.contains("NH2") || ms2Ion.ionLabel.contains("H2O") || ms2Ion.ionLabel.contains("a");
//        };
//
//        const auto termnator = std::remove_if(ms2Ions->begin(), ms2Ions->end(), terminatorLogic);
//        ms2Ions->erase(termnator, ms2Ions->end());
//    }
//
//}//namespace
//void FragLibReader::getTopNMostIntenseMs2Ions(
//        int topNMs2Ions,
//        double mzMin,
//        double mzMax,
//        bool byIonsOnly,
//        QVector<MS2Ion> *ms2Ions
//) {
//
//    removeMzVals(mzMin, mzMax, ms2Ions);
//
//    if (byIonsOnly) {
//        removeNeutralLosses(ms2Ions);
//    }
//
//    const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
//        return l.intensity < r.intensity;
//    };
//
//    std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);
//
//    topNMs2Ions = std::min(topNMs2Ions, ms2Ions->size());
//
//    ms2Ions->resize(topNMs2Ions);
//
//    for (int i = 0; i < ms2Ions->size(); i++) {
//        (*ms2Ions)[i].rank = i;
//    }
//
//}
//
//void FragLibReader::filterMs2IonsByIntensity(
//        double intensityThreshold,
//        QVector<MS2Ion> *ms2Ions
//        ) {
//
//    const auto terminatorLogic = [intensityThreshold](const MS2Ion &ion){
//        return ion.intensity < intensityThreshold;
//    };
//
//    const auto terminator = std::remove_if(
//            ms2Ions->begin(),
//            ms2Ions->end(),
//            terminatorLogic
//            );
//
//    ms2Ions->erase(terminator, ms2Ions->end());
//
//}
//
//namespace {
//
//    double calculateMassFromThomson(
//            double mz,
//            int charge
//    ) {
//        const double PROTON = 1.0072;
//        return (mz * charge) - (charge * PROTON);
//    }
//
//}//namespace
//Err FragLibReader::buildFragIonLibForCandidates(
//        const QString &fragLibUri,
//        int chargeMin,
//        int chargeMax,
//        double mzTargetMin,
//        double mzTargetMax,
//        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
//        ) {
//
//    ERR_INIT
//
////    FragLibReader fragLibReader;
////    e = fragLibReader.init(fragLibUri); ree;
////
////    const int topNMs2Ions = 1000;
////
////    for (Charge charge = chargeMin; charge <= chargeMax; ++charge) {
////
////        const double massStart = calculateMassFromThomson(
////                mzTargetMin,
////                charge
////        );
////
////        const double massEnd = calculateMassFromThomson(
////                mzTargetMax,
////                charge
////        );
////
////        const double mzMin = 100.0;
////        const double mzMax = 2000.0;
////
////        e = fragLibReader.getMS2Ions(
////                massStart,
////                massEnd,
////                mzMin,
////                mzMax,
////                topNMs2Ions,
////                peptideSequenceChargeKeyVsCandidatePeptide
////        ); ree
//
////    }
//
//    ERR_RETURN
//}
//






