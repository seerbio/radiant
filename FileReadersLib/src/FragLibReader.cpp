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
        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideSequenceChargeKeyVsMS2Ions,
        QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsMass,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsIRT
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
            peptideSequenceChargeKeyVsIsDecoy,
            peptideSequenceChargeKeyVsMass,
            peptideSequenceChargeKeyVsIRT
            ); ree;

    qDebug() << "MS2 Predictions count:" << peptideSequenceChargeKeyVsMS2Ions->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN
}

Err FragLibReader::getMS2Ions(
        double massStart,
        double massEnd,
        QMap<PeptideSequenceChargeKey, MS2IonsSeparated> *peptideSequenceChargeKeyVsMS2IonsSeparated,
        QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsMass,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsIRT
) {

    ERR_INIT

    QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> peptideSequenceChargeKeyVsMS2Ions;
    e = getMS2Ions(
            massStart,
            massEnd,
            &peptideSequenceChargeKeyVsMS2Ions,
            peptideSequenceChargeKeyVsIsDecoy,
            peptideSequenceChargeKeyVsMass,
            peptideSequenceChargeKeyVsIRT
            ); ree

    e = buildMS2IonsSeparated(
            peptideSequenceChargeKeyVsMS2Ions,
            peptideSequenceChargeKeyVsMS2IonsSeparated
            ); ree;

    ERR_RETURN
}

Err FragLibReader::getMS2Ions(
        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideSequenceChargeKeyVsMS2Ions,
        QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsMass,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsIRT
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
            peptideSequenceChargeKeyVsMS2Ions,
            peptideSequenceChargeKeyVsIsDecoy,
            peptideSequenceChargeKeyVsMass,
            peptideSequenceChargeKeyVsIRT
    ); ree

    qDebug() << "All MS2 Predictions count:" << peptideSequenceChargeKeyVsMS2Ions->size() << "retrieved in" << et.elapsed() << "mSec";

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
        QMap<PeptideStringWithMods, MS2IonsSeparated> *fragPreds,
        QMap<PeptideStringWithMods, bool> *fragPredsIsDecoy,
        QMap<PeptideStringWithMods, double> *fragPredsMass,
        QMap<PeptideStringWithMods, double> *fragPredsIRT
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

        QMap<PeptideSequenceChargeKey, MS2IonsSeparated> peptideSequenceChargeKeyVsMS2IonsSeparated;
        QMap<PeptideSequenceChargeKey, bool> peptideSequenceChargeKeyVsIsDecoy;
        QMap<PeptideSequenceChargeKey, double> peptideSequenceChargeKeyVsMass;
        QMap<PeptideSequenceChargeKey, double> peptideSequenceChargeKeyVsIRT;

        e = fragLibReader.getMS2Ions(
                massStart,
                massEnd,
                &peptideSequenceChargeKeyVsMS2IonsSeparated,
                &peptideSequenceChargeKeyVsIsDecoy,
                &peptideSequenceChargeKeyVsMass,
                &peptideSequenceChargeKeyVsIRT
        ); ree

        for (
            auto it = peptideSequenceChargeKeyVsMS2IonsSeparated.begin();
            it != peptideSequenceChargeKeyVsMS2IonsSeparated.end();
            it++
            ) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            const MS2IonsSeparated &ms2IonsSeparated = it.value();

            PeptideStringWithMods peptideStringWithMods;
            Charge chargeFromPeptideSequenceChargeKey;
            e = peptideStringWithModsFromPeptideSequenceChargeKey(
                    peptideSequenceChargeKey,
                    &peptideStringWithMods,
                    &chargeFromPeptideSequenceChargeKey
            ); ree;

            if (charge != chargeFromPeptideSequenceChargeKey) {
                continue;
            }

            e = ErrorUtils::isTrue(peptideSequenceChargeKeyVsIsDecoy.contains(peptideSequenceChargeKey)); ree;

            fragPreds->insert(peptideStringWithMods, ms2IonsSeparated);
            fragPredsIsDecoy->insert(peptideStringWithMods, peptideSequenceChargeKeyVsIsDecoy.value(peptideSequenceChargeKey));
            fragPredsMass->insert(peptideStringWithMods, peptideSequenceChargeKeyVsMass.value(peptideSequenceChargeKey));
            fragPredsIRT->insert(peptideStringWithMods, peptideSequenceChargeKeyVsIRT.value(peptideSequenceChargeKey));
        }

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
        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> *peptideStringWithModsVsMS2Ions,
        QMap<PeptideSequenceChargeKey, bool> *peptideSequenceChargeKeyVsIsDecoy,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsMass,
        QMap<PeptideSequenceChargeKey, double> *peptideSequenceChargeKeyVsIRT
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

        peptideSequenceChargeKeyVsMass->insert(
                row.peptideSequenceChargeKey,
                row.mass
        );

        peptideSequenceChargeKeyVsIRT->insert(
                row.peptideSequenceChargeKey,
                row.iRT
        );
    }

    ERR_RETURN
}

Err FragLibReader::buildMS2IonsSeparated(
        const QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> &peptideSequenceChargeKeyVsMS2Ions,
        QMap<PeptideSequenceChargeKey, MS2IonsSeparated> *peptideSequenceChargeKeyVsMS2IonsSeparated
        ) {

    ERR_INIT

    for (auto it = peptideSequenceChargeKeyVsMS2Ions.begin();
         it != peptideSequenceChargeKeyVsMS2Ions.end(); it++) {

        const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
        const QVector<MS2Ion> &ms2Ions = it.value();

        MS2IonsSeparated ms2IonsSeparated;
        for (const MS2Ion &ion : ms2Ions) {

            if (ion.mz < 0) {
                continue; //TODO make sure no sub zeros are entered when building the library in IstrosLibBuilder workflow
            }

            QPair<IonIndex, IonType> ionInfo;
            e = ion.getIonLabelInfo(&ionInfo); ree

            if (ionInfo.second == "y") {
                ms2IonsSeparated.yIons.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "y-H2O") {
                ms2IonsSeparated.yH2OIons.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "y-NH3") {
                ms2IonsSeparated.yNH3Ions.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "y^2") {
                ms2IonsSeparated.y2Ions.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "b") {
                ms2IonsSeparated.bIons.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "b-H2O") {
                ms2IonsSeparated.bH2OIons.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "b-NH3") {
                ms2IonsSeparated.bNH3Ions.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "b^2") {
                ms2IonsSeparated.b2Ions.insert(ionInfo.first, ion);
            }
            else if (ionInfo.second == "a") {
                ms2IonsSeparated.aIons.insert(ionInfo.first, ion);
            }
            else {
                rrr(eValueError);
            }
        }

        peptideSequenceChargeKeyVsMS2IonsSeparated->insert(peptideSequenceChargeKey, ms2IonsSeparated);

    }

    ERR_RETURN
}

namespace {

    RTree loadRtree(const QMap<PeptideString, QVector<MS2Ion>> &peptideSequenceVsMS2Ions) {

        std::vector<rTreePoint> cloudLoader;

        for (auto it = peptideSequenceVsMS2Ions.begin(); it != peptideSequenceVsMS2Ions.end(); it++) {

            const PeptideString &peptideStr = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value();

            for (const MS2Ion &ms2Ion : ms2Ions) {
                rTreeCoor coor(ms2Ion.mz, 0.0);
                cloudLoader.emplace_back(coor, -1);
            }
        }

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    QMap<MzHashed, MZION> buildMzHashedVsMzIon(const QList<QVector<MS2Ion>> &ms2IonVecs) {

        QMap<MzHashed, MZION> mzHashedVsMzIon;

        const int precision = S_GLOBAL_SETTINGS.HASHING_PRECISION;

        for (const QVector<MS2Ion> &ms2Ions : ms2IonVecs) {

            for (const MS2Ion &ms2Ion : ms2Ions) {

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, precision);
                const MZION mzIon = MathUtils::pRound(ms2Ion.mz, precision);

                mzHashedVsMzIon.insert(mzHashed, mzIon);
            }
        }

        return mzHashedVsMzIon;
    }

}
Err FragLibReader::generateFragmentFrequencies(
        const QMap<PeptideString, QVector<MS2Ion>> &peptideStringVsMS2Ions,
        double ppmTol,
        QMap<MzHashed, FrequencyPercent> *fragmentFrequencies
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(peptideStringVsMS2Ions); ree
    fragmentFrequencies->clear();

    RTree rTree = loadRtree(peptideStringVsMS2Ions);
    const QMap<MzHashed, MZION> mzIonHashedVsMzIon = buildMzHashedVsMzIon(peptideStringVsMS2Ions.values());

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

Err FragLibReader::unseparateMS2IonsSeparated(
        const QMap<PeptideString, MS2IonsSeparated> &peptideSequenceVsMS2IonsSeparated,
        QMap<PeptideString, QVector<MS2Ion>> *peptideSequenceVsMS2Ions) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(peptideSequenceVsMS2IonsSeparated); ree

    for (auto it = peptideSequenceVsMS2IonsSeparated.begin(); it != peptideSequenceVsMS2IonsSeparated.end(); it++) {

        const PeptideString &peptideString = it.key();
        const MS2IonsSeparated &ms2IonsSeparated = it.value();

        QVector<MS2Ion> ms2IonsConcat;
        ms2IonsConcat.append(ms2IonsSeparated.yIons.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.y2Ions.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.bIons.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.b2Ions.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.yNH3Ions.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.yH2OIons.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.bNH3Ions.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.bH2OIons.values().toVector());
        ms2IonsConcat.append(ms2IonsSeparated.aIons.values().toVector());

        peptideSequenceVsMS2Ions->insert(peptideString, ms2IonsConcat);
    }

    ERR_RETURN
}
