//
// Created by anichols on 4/3/23.
//

#include "FragLibReader.h"

#include "LibraryCommon.h"
#include "LibraryReader.h"
#include "MolecularFormula.h"
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
, m_mzMax(2000.0)
, m_useBYOnly(false)
{}

Err FragLibReader::init(const QString &fragLibFilePath, const AminoAcids &aminoAcids) {
    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    m_fragLibFilePath = fragLibFilePath;

    e = ParquetReader::read(
            m_fragLibFilePath,
            &m_fragLibReaderRows
    ); ree

    m_aminoAcids = aminoAcids;

    qDebug() << "Use b/y ions only" << m_useBYOnly;

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
            ion.charge = ion.ionLabel.contains("^2") ? 2 : 1;
            ms2Ions->push_back(ion);
        }

        ERR_RETURN
    }

}//namespace
Err FragLibReader::getMS2Ions(
        const QString &fragLibFilePath,
        const AminoAcids &aminoAcids,
        double massStart,
        double massEnd,
        double mzMin,
        double mzMax,
        int topNMs2Ions,
        bool byIonsOnly,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = ParquetReader::read(
            fragLibFilePath,
            FragLibReaderNamespace::MASS,
            {massStart, massEnd},
            &fragLibReaderRows
            ); ree;

    e = fragLibReaderRowsToMs2IonsMap(
            fragLibReaderRows,
            aminoAcids,
            topNMs2Ions,
            mzMin,
            mzMax,
            byIonsOnly,
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
            m_mzMin,
            m_mzMax,
            m_useBYOnly,
            peptideSequenceChargeKeyVsCandidatePeptide
            ); ree;

    ERR_RETURN
}

Err FragLibReader::getMS2IonsTopN(
        int topNMs2Ions,
        double mzMin,
        double mzMax,
        bool byIonsOnly,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = fragLibReaderRowsToMs2IonsMap(
            m_fragLibReaderRows,
            m_aminoAcids,
            topNMs2Ions,
            mzMin,
            mzMax,
            byIonsOnly,
            peptideSequenceChargeKeyVsCandidatePeptide
    ); ree

    qDebug() << "All MS2 Predictions count:" << peptideSequenceChargeKeyVsCandidatePeptide->size() << "retrieved in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err FragLibReader::getMS2IonsTopN(
        const QMap<Index, bool> &selectionList,
        int topNMs2Ions,
        double mzMin,
        double mzMax,
        bool byIonsOnly,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
        ) {

    ERR_INIT

    QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptideTemp;
    e = getMS2IonsTopN(
            topNMs2Ions,
            mzMin,
            mzMax,
            byIonsOnly,
            &peptideSequenceChargeKeyVsCandidatePeptideTemp
            ); ree;

    int counter = 0;
    for (auto it = peptideSequenceChargeKeyVsCandidatePeptideTemp.begin(); it != peptideSequenceChargeKeyVsCandidatePeptideTemp.end(); it++) {

        const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
        const CandidatePeptide &candidatePeptide = it.value();

        //for libraries w/ decoys included.
//        if (candidatePeptide.isDecoy) {
//            continue;
//        }

        if (selectionList.value(counter++)) {
            peptideSequenceChargeKeyVsCandidatePeptide->insert(peptideSequenceChargeKey, candidatePeptide);
        }

    }

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

    void removeNeutralLosses(QVector<MS2Ion> *ms2Ions) {

        const auto terminatorLogic = [](const MS2Ion &ms2Ion){
            return ms2Ion.ionLabel.contains("NH2") || ms2Ion.ionLabel.contains("H2O") || ms2Ion.ionLabel.contains("a");
        };

        const auto termnator = std::remove_if(ms2Ions->begin(), ms2Ions->end(), terminatorLogic);
        ms2Ions->erase(termnator, ms2Ions->end());
    }

}//namespace
void FragLibReader::getTopNMostIntenseMs2Ions(
        int topNMs2Ions,
        double mzMin,
        double mzMax,
        bool byIonsOnly,
        QVector<MS2Ion> *ms2Ions
) {

    removeMzVals(mzMin, mzMax, ms2Ions);

    if (byIonsOnly) {
        removeNeutralLosses(ms2Ions);
    }

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

//    FragLibReader fragLibReader;
//    e = fragLibReader.init(fragLibUri); ree;
//
//    const int topNMs2Ions = 1000;
//
//    for (Charge charge = chargeMin; charge <= chargeMax; ++charge) {
//
//        const double massStart = calculateMassFromThomson(
//                mzTargetMin,
//                charge
//        );
//
//        const double massEnd = calculateMassFromThomson(
//                mzTargetMax,
//                charge
//        );
//
//        const double mzMin = 100.0;
//        const double mzMax = 2000.0;
//
//        e = fragLibReader.getMS2Ions(
//                massStart,
//                massEnd,
//                mzMin,
//                mzMax,
//                topNMs2Ions,
//                peptideSequenceChargeKeyVsCandidatePeptide
//        ); ree

//    }

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

namespace {

    double calculateFragmentMass(
            const QString &fragment,
            const AminoAcids &aminoAcids,
            bool isCTerm
            ){

        const auto accLogic = [aminoAcids](double sum, const QChar &c){
            return sum + aminoAcids.aminoAcid(c).monoisotopicMass();
        };

        const double fragMass = std::accumulate(fragment.begin(), fragment.end(), 0.0, accLogic);
        return isCTerm ? fragMass + Molecule(MolecularFormulas::waterFormula).monoisotopicMass() : fragMass;
    }

}//namespace
Err FragLibReader::fragLibReaderRowsToMs2IonsMap(
        const QVector<FragLibReaderRow> &fragLibReaderRows,
        const AminoAcids &aminoAcids,
        int topNMs2Ions,
        double mzMin,
        double mzMax,
        bool byIonsOnly,
        QMap<PeptideSequenceChargeKey, CandidatePeptide> *peptideSequenceChargeKeyVsCandidatePeptide
) {

    ERR_INIT

    for (const FragLibReaderRow &row : fragLibReaderRows) {

//        e = ErrorUtils::isFalse(peptideSequenceChargeKeyVsCandidatePeptide->contains(row.peptideSequenceChargeKey));
//        if (e != eNoError) {
////            qDebug() << row.peptideSequenceChargeKey; //TODO correct this.
////            rrr(eValueError);
//        }

        PeptideStringWithMods peptideStringWithMods;
        Charge charge;
        e = FragLibReader::peptideStringWithModsFromPeptideSequenceChargeKey(
                row.peptideSequenceChargeKey,
                &peptideStringWithMods,
                &charge
                ); ree;

        QVector<MS2Ion> ms2Ions;
        e = buildMzIons(
                row,
                &ms2Ions
        ); ree;

        CandidatePeptide candidatePeptide;
        candidatePeptide.totalFragmentCount = ms2Ions.size();
        candidatePeptide.peptideStringWithMods = peptideStringWithMods;

        getTopNMostIntenseMs2Ions(
                topNMs2Ions,
                mzMin,
                mzMax,
                byIonsOnly,
                &ms2Ions
                ); ree

        candidatePeptide.ms2Ions = ms2Ions;

        candidatePeptide.ms2IonMzB2B3 = {
                calculateFragmentMass(peptideStringWithMods.mid(0, 2), aminoAcids, false),
                calculateFragmentMass(peptideStringWithMods.mid(0, 3), aminoAcids, false)
        };

        candidatePeptide.mass = row.mass;
        candidatePeptide.iRt = row.iRT;
        candidatePeptide.isDecoy = row.isDecoy;
        candidatePeptide.charge = charge;

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

//    qDebug() << "Frequences calculated in" << et.elapsed() << "mSec";

    ERR_RETURN
}


namespace {

    //From biophysicalcalcs.h
    double calculatePeptideMass(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            const QHash<ResidueIndex, ModificationMass> &mods
    ) {

        const auto accumulateLogic
                = [&](double sum, const QChar &aa){return sum + aminoAcids.aminoAcid(aa).monoisotopicMass(); };

        const double mass = std::accumulate(sequence.begin(),  sequence.end(), CommonMolecules::H2O.monoisotopicMass(), accumulateLogic);

        const QList<double> &modMasses = mods.values();
        const double modMassesSum = std::accumulate(modMasses.begin(), modMasses.end(), 0.0);

        return mass + modMassesSum;
    }

    void removeModificationLabel(PeptideStringWithMods *peptideStringWithMods) {

        bool bypassOn = false;

        PeptideStringWithMods peptideStringWithModsStripped;
        for (const QChar &aa : *peptideStringWithMods) {

            if (aa == '(') {
                bypassOn = true;
            }
            else if (aa == ')') {
                bypassOn = false;
                continue;
            }

            if (bypassOn) {
                continue;
            }

            peptideStringWithModsStripped += aa;
        }

        *peptideStringWithMods = peptideStringWithModsStripped;
    }


}//namespace
Err FragLibReader::convertDIANNLibToFragLib(
        const QString &specLibFilePath,
        const AminoAcids &aminoAcids
        ) {

    ERR_INIT

//    AminoAcids aminoAcids;
//    aminoAcids.addFixedModification('C', MolecularFormulas::carbamidomethylFormula);

    LibraryReader lib;
    lib.load(specLibFilePath.toStdString());

    QVector<FragLibReaderRow> fragLibReaderRows;

    for (const Entry &entry : lib.getEntries()) {

        const QString diannPepSeqChargeString = QString::fromStdString(entry.getName());

        Charge charge;
        e = ErrorUtils::toInt(QString(diannPepSeqChargeString.back()), &charge); ree;

        //TODO reexplore this. So you want to use charge 4's?
        if (charge == 4 || charge == 1) {
            continue;
        }

        PeptideStringWithMods peptideStringWithMods = diannPepSeqChargeString;
        peptideStringWithMods = peptideStringWithMods.replace(QString::number(charge), "");
        removeModificationLabel(&peptideStringWithMods);

        const PeptideSequenceChargeKey peptideSequenceChargeKey = peptideStringWithMods + "|" + QString::number(charge);
        const double mass = calculatePeptideMass(peptideStringWithMods, aminoAcids, {});
        const double iRT = entry.getTarget().getIRT();

        QVector<double> mzVals;
        QVector<double> intensityVals;
        QStringList ionLabels;
        for (const Product &pr : entry.getTarget().getFragments()) {

            const int mzCharge = toascii(static_cast<unsigned char>(pr.charge));
            const QChar ionType = toascii(pr.type) == 1 ? 'b' : 'y';
            const int index = ionType == 'b' ? toascii(pr.index) + 1 : peptideStringWithMods.size() - toascii(pr.index);
            const QString ionModifier = mzCharge == 2
                    ?  "^2" : QString::fromStdString(LibraryCommon::getLoss(static_cast<int>(toascii(pr.loss))));

            const QString ionLabel = ionType + QString::number(index) + ionModifier;
            const double mz = pr.mz;
            const double intensity = pr.height;

            mzVals.push_back(mz);
            intensityVals.push_back(intensity);
            ionLabels.push_back(ionLabel);
        }

        const QString ionLabelsJoined = ionLabels.join(";");

        FragLibReaderRow row;
        row.peptideSequenceChargeKey = peptideSequenceChargeKey;
        row.mzVals = mzVals;
        row.intensityVals = intensityVals;
        row.ionLabels = ionLabelsJoined;
        row.mass = mass;
        row.isDecoy = 0;
        row.iRT = iRT;

        fragLibReaderRows.push_back(row);
    }

    const QString outputFilePath = specLibFilePath + ".fragLib";
    e = ParquetReader::write(fragLibReaderRows, outputFilePath); ree;

    ERR_RETURN
}

int FragLibReader::libarySize() {
    return m_fragLibReaderRows.size();
}

Err FragLibReader::mutateCandidatePeptideTarget(
        const CandidatePeptide &candidatePeptideTarget,
        CandidatePeptide *candidatePeptideDecoy
        ) {

    ERR_INIT

    const QMap<QChar, double> diannMutateAminoAcidTo = AminoAcids::diannMutateAminoAcidTo();

    const QString &seq = candidatePeptideTarget.peptideStringWithMods;

    const int firstIndexToMutate = 1;
    const int secondIndexToMutate = seq.size() - 2;

    const double nTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(firstIndexToMutate));
    const double cTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(secondIndexToMutate));
    const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
    const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

    QVector<MS2Ion> ms2IonDecoys;
    for (const MS2Ion &ms2Ion : candidatePeptideTarget.ms2Ions) {

        MS2Ion ms2IonDecoy = ms2Ion;

        QPair<IonIndex, IonType> ionLableInfo;
        e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

        if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

            if (ionLableInfo.second.contains("^2")) {

                if (ionLableInfo.first >= firstIndexToMutate) {
                    ms2IonDecoy.mz += nTermDeltaMassCharge2;
                }

                if (ionLableInfo.first >= secondIndexToMutate) {
                    ms2IonDecoy.mz += cTermDeltaMassCharge2;
                }
            }
            else {

                if (ionLableInfo.first >= firstIndexToMutate) {
                    ms2IonDecoy.mz += nTermDeltaMass;
                }

                if (ionLableInfo.first >= secondIndexToMutate) {
                    ms2IonDecoy.mz += cTermDeltaMass;
                }
            }
        }

        else if (ionLableInfo.second.contains('y')) {

            if (ionLableInfo.second.contains("^2")) {

                if (ionLableInfo.first >= firstIndexToMutate) {
                    ms2IonDecoy.mz += cTermDeltaMassCharge2;
                }

                if (ionLableInfo.first >= secondIndexToMutate) {
                    ms2IonDecoy.mz += nTermDeltaMassCharge2;
                }
            }
            else {

                if (ionLableInfo.first >= firstIndexToMutate) {
                    ms2IonDecoy.mz += cTermDeltaMass;
                }

                if (ionLableInfo.first >= secondIndexToMutate) {
                    ms2IonDecoy.mz += nTermDeltaMass;
                }
            }
        }

        else {
            qDebug() << "Non b/y/a ion" << ionLableInfo;
        }

        ms2IonDecoys.push_back(ms2IonDecoy);
    }

    QVector<MZION> ms2IonMzB2B3Decoy;
    std::transform(
            candidatePeptideTarget.ms2IonMzB2B3.begin(),
            candidatePeptideTarget.ms2IonMzB2B3.end(),
            std::back_inserter(ms2IonMzB2B3Decoy),
            [nTermDeltaMass](MZION mz){return mz + nTermDeltaMass;}
            );

    *candidatePeptideDecoy = candidatePeptideTarget;
    candidatePeptideDecoy->isDecoy = true;
    candidatePeptideDecoy->ms2Ions = ms2IonDecoys;
    candidatePeptideDecoy->ms2IonMzB2B3 = ms2IonMzB2B3Decoy;

    ERR_RETURN
}


