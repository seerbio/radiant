//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePairManager.h"

#include "AminoAcids.h"
#include "BiophysicalCalcs.h"
#include "FragLibReader.h"
#include "MolecularFormula.h"
#include "Molecule.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class Q_DECL_HIDDEN TargetDecoyCandidatePairManager::Private
{
    using RTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using RTreeBox = bg::model::box<RTreeCoor>;
    using RTreePoint = std::pair<RTreeCoor, TargetDecoyCandidatePair*> ;
    using RTree = bgi::rtree<RTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(QVector<TargetDecoyCandidatePair>* targetDecoyCandidatePairPntrs);

    Err getPoints(
            double mzMin,
            double mzMax,
            QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
    );

private:

    bool m_isInit;
    RTree *m_rTree;

};

TargetDecoyCandidatePairManager::Private::Private()
        : m_rTree(Q_NULLPTR)
        , m_isInit(false)
{}

TargetDecoyCandidatePairManager::Private::~Private() {
    delete m_rTree;
}

Err TargetDecoyCandidatePairManager::Private::init(QVector<TargetDecoyCandidatePair> *targetDecoyCandidatePairPntrs) {

    ERR_INIT

    e = ErrorUtils::isFalse(targetDecoyCandidatePairPntrs->isEmpty()); ree;

    std::vector<RTreePoint> cloudLoader;
    cloudLoader.reserve(targetDecoyCandidatePairPntrs->size());
    std::transform(
            targetDecoyCandidatePairPntrs->begin(),
            targetDecoyCandidatePairPntrs->end(),
            std::back_inserter(cloudLoader),
            [](TargetDecoyCandidatePair &tdcp){
                const RTreeCoor coor(tdcp.mz(), 0.0);
                return RTreePoint(coor, &tdcp);
            }
    );

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    m_isInit = true;
    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::Private::getPoints(
        double mzMin,
        double mzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
) {

    ERR_INIT

    e = ErrorUtils::isTrue(mzMin <= mzMax); ree;
    e = ErrorUtils::isTrue(m_isInit); ree;

    const RTreeBox queryBox(
            RTreeCoor(mzMin, 0.0),
            RTreeCoor(mzMax, 0.0)
    );

    std::vector<RTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const RTreePoint &rtp : result) {
        targetDecoyCandidatePairs->push_back(rtp.second);
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


TargetDecoyCandidatePairManager::TargetDecoyCandidatePairManager()
: d_ptr(new Private())
{}

TargetDecoyCandidatePairManager::~TargetDecoyCandidatePairManager() {}

namespace {

    void removeNonSequencesWithNonCanonicalAminoAcids(QVector<FragLibReaderRow> *fragLibReaderRows) {

        const QString verbodenAminoAcids = QStringLiteral("BJOUXZ");
        const auto terminatorLogic = [verbodenAminoAcids](const FragLibReaderRow &flrr){
            return std::any_of(
                    verbodenAminoAcids.begin(),
                    verbodenAminoAcids.end(),
                    [flrr](const QChar &aa){return flrr.peptideSequenceChargeKey.contains(aa);}
                    ) ;
        };

        const auto terminator = std::remove_if(
                fragLibReaderRows->begin(),
                fragLibReaderRows->end(),
                terminatorLogic
                );

        fragLibReaderRows->erase(terminator, fragLibReaderRows->end());
    }

}
Err TargetDecoyCandidatePairManager::init(
        const PythiaParameters &pythiaParameters,
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isFalse(fragLibReaderRows->isEmpty()); ree;

    m_pythiaParameters = pythiaParameters;

    removeNonSequencesWithNonCanonicalAminoAcids(fragLibReaderRows);

    e = buildTargetDecoyCandidatePairs(*fragLibReaderRows); ree;

    e = d_ptr->init(&m_targetDecoyCandidatePairs); ree;

    ERR_RETURN
}

namespace {

    Err buildMS2Ions(
            const PythiaParameters &pythiaParameters,
            const FragLibReaderRow &flrr,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(flrr.mzVals); ree;
        e = ErrorUtils::isEqual(flrr.mzVals.size(), flrr.intensityVals.size());ree;
        e = ErrorUtils::isNotEmpty(flrr.ionLabels); ree

        const QStringList ionLabelsSplit = flrr.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, Qt::SkipEmptyParts);
        e = ErrorUtils::isEqual(flrr.mzVals.size(), ionLabelsSplit.size());ree;

        QVector<MS2Ion> ms2IonsBuilder;
        ms2IonsBuilder.reserve(flrr.mzVals.size());

        for (int i = 0; i < flrr.mzVals.size(); i++) {
            MS2Ion ms2Ion;
            ms2Ion.mz = flrr.mzVals.at(i);
            ms2Ion.intensity = flrr.intensityVals.at(i);
            ms2Ion.ionLabel = ionLabelsSplit.at(i);
            ms2Ion.iRT = static_cast<IRT>(flrr.iRT);
            ms2Ion.charge = ms2Ion.ionLabel.contains("^2") ? 2 : 1;

            ms2IonsBuilder.push_back(ms2Ion);
        }

        MS2Ion::filterMS2IonsByMz(
                pythiaParameters.mzMinDataStructure,
                pythiaParameters.mzMaxDataStructure,
                &ms2IonsBuilder
        );

        MS2Ion::sortMS2IonsIntensityDesc(&ms2IonsBuilder);

        const int ms2IonsSize = std::min(pythiaParameters.topNMs2Ions, ms2IonsBuilder.size());
        ms2Ions->reserve(ms2IonsSize);
        ms2IonsBuilder.resize(ms2IonsSize);

        *ms2Ions = ms2IonsBuilder;

        ERR_RETURN
    }

    Err mutateCandidatePeptideTarget(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonTarget,
            QVector<MS2Ion> *ms2IonDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;
        e = ErrorUtils::isNotEmpty(ms2IonTarget); ree;

        ms2IonDecoys->clear();
        ms2IonDecoys->reserve(ms2IonTarget.size());

        const QMap<QChar, double> diannMutateAminoAcidTo = AminoAcids::diannMutateAminoAcidTo();

        const QString &seq = peptideStringWithMods;

        const int firstIndexToMutate = 1;
        const int secondIndexToMutate = seq.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidTo.value(seq.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

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

            ms2IonDecoys->push_back(ms2IonDecoy);
        }

        ERR_RETURN
    }

    Err reverseCandidatePeptideTarget(
            const PeptideStringWithMods &peptideStringWithMods,
            const AminoAcids &aminoAcids,
            const QVector<MS2Ion> &ms2IonTarget,
            QVector<MS2Ion> *ms2IonDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;
        e = ErrorUtils::isNotEmpty(ms2IonTarget); ree;

        ms2IonDecoys->clear();
        ms2IonDecoys->reserve(ms2IonTarget.size());

        PeptideStringWithMods peptideStringWithModsMiddelReversed = peptideStringWithMods;
        std::reverse(
                peptideStringWithModsMiddelReversed.begin() + 1,
                peptideStringWithModsMiddelReversed.begin() + peptideStringWithModsMiddelReversed.size() - 1
                );

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo); ree;

            if (ionLableInfo.second.contains('b')) {
                const QString bSeq = peptideStringWithModsMiddelReversed.left(ionLableInfo.first);
                ms2IonDecoy.mz = BiophysicalCalcs::calculateThomson(bSeq, aminoAcids, ms2IonDecoy.charge);
            }

            else if (ionLableInfo.second.contains('y')) {

                const QString ySeq = peptideStringWithModsMiddelReversed.right(ionLableInfo.first);
                ms2IonDecoy.mz
                    = BiophysicalCalcs::calculateThomson(ySeq, aminoAcids, ms2IonDecoy.charge) +  (CommonMolecules::H2O.monoisotopicMass() / ms2Ion.charge);
            }

            else {
                qDebug() << "Non b/y/a ion" << ionLableInfo;
            }

            ms2IonDecoys->push_back(ms2IonDecoy);
        }

        ERR_RETURN
    }

    QPair<Err, TargetDecoyCandidatePair> targetDecoyCandidatePairsLoadLogic(
            const FragLibReaderRow &flrr,
            const PythiaParameters &pythiaParameters
            ) {

        ERR_INIT

        PeptideStringWithMods peptideStringWithMods;
        Charge charge;
        e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
                flrr.peptideSequenceChargeKey,
                &peptideStringWithMods,
                &charge
        ); rree;

        //TODO once incorporating PTMs in the sequence, make algo to remove them
        // and then calculate the peptide length MODS
        const int peptideLength = peptideStringWithMods.size();
        if (peptideLength < pythiaParameters.peptideLengthMin
            || peptideLength > pythiaParameters.peptideLengthMax
            || charge < pythiaParameters.chargeStateMin
            || charge > pythiaParameters.chargeStateMax
            || flrr.isDecoy) {
            return {e, {}};
        }

        QVector<MS2Ion> ms2IonsTarget;
        e = buildMS2Ions(
                pythiaParameters,
                flrr,
                &ms2IonsTarget
        ); rree;

        QVector<MS2Ion> ms2IonsDecoy;

//#define REVERSE_MIDDLE_DECOY
#ifdef REVERSE_MIDDLE_DECOY
        e= reverseCandidatePeptideTarget(
                peptideStringWithMods,
                pythiaParameters.aminoAcids,
                ms2IonsTarget,
                &ms2IonsDecoy
        ); rree;
#else
        e= mutateCandidatePeptideTarget(
                peptideStringWithMods,
                ms2IonsTarget,
                &ms2IonsDecoy
        ); rree;
#endif

        TargetDecoyCandidatePair targetDecoyCandidatePair(
                peptideStringWithMods,
                ms2IonsTarget,
                ms2IonsDecoy,
                flrr.precursorCharge,
                flrr.mass,
                flrr.iRT,
                flrr.mzVals.size()
        );

        return {e, targetDecoyCandidatePair};
    }

}//namespace
Err TargetDecoyCandidatePairManager::buildTargetDecoyCandidatePairs(
        const QVector<FragLibReaderRow> &fragLibReaderRows
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibReaderRows); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;

    QElapsedTimer et;
    et.start();

    m_targetDecoyCandidatePairs.reserve(fragLibReaderRows.size()); ree;

#define PARALLEL_FRAGLIB_LOAD
#ifdef PARALLEL_FRAGLIB_LOAD
    const auto loadLogicBinder = std::bind(
            targetDecoyCandidatePairsLoadLogic,
            std::placeholders::_1,
            m_pythiaParameters
    );

    QFuture<QPair<Err, TargetDecoyCandidatePair>> futures = QtConcurrent::mapped(
            fragLibReaderRows,
            loadLogicBinder
            );
    futures.waitForFinished();

    for (const QPair<Err, TargetDecoyCandidatePair> &result : futures) {
        e = result.first; ree;
        const TargetDecoyCandidatePair &tdcp = result.second;
        m_targetDecoyCandidatePairs.push_back(tdcp);
    }
#else
    for (const FragLibReaderRow &flrr : fragLibReaderRows) {

        QPair<Err, TargetDecoyCandidatePair> result = targetDecoyCandidatePairsLoadLogic(
                flrr,
                m_pythiaParameters
                );

        e = result.first; ree;

        const TargetDecoyCandidatePair &tdcp = result.second;
        m_targetDecoyCandidatePairs.push_back(tdcp);
    }
#endif

    qDebug() << m_targetDecoyCandidatePairs.size() << "Candidates loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(
        double mzMin,
        double mzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isTrue(mzMax >= mzMin); ree;

    e = d_ptr->getPoints(
            mzMin,
            mzMax,
            targetDecoyPointers
            ); ree;

    std::mt19937 gen(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
    std::shuffle(targetDecoyPointers->begin(), targetDecoyPointers->end(), gen);

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::getTargetDecoyCandidatePairPointers(
        double mzMin,
        double mzMax,
        double randomSelectionFraction,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(randomSelectionFraction <= 1.0); ree;

    QVector<TargetDecoyCandidatePair*> targetDecoyPointersAll;
    e = getTargetDecoyCandidatePairPointers(
            mzMin,
            mzMax,
            &targetDecoyPointersAll
            ); ree;

    if (randomSelectionFraction < 0) {
        *targetDecoyPointers = targetDecoyPointersAll;
        ERR_RETURN
    }

    const int testDataSize = static_cast<int>(targetDecoyPointersAll.size() * randomSelectionFraction);
    const int seed = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;

    const QMap<int, bool> selectionList = MathUtils::generateRandomSelectionList(
            targetDecoyPointersAll.size(),
            testDataSize,
            seed
            );

    for (int i = 0; i < targetDecoyPointersAll.size(); i++) {

        if (!selectionList.value(i)) {
            continue;
        }

        targetDecoyPointers->push_back(targetDecoyPointersAll[i]);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        PeptideStringWithMods *peptideStringWithMods,
        Charge *charge
){

    ERR_INIT

    const int expectedSplitSize = 2;

    const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
            S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
            Qt::SkipEmptyParts
    );

    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize
            ); ree;

    *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
    ); ree

    ERR_RETURN
}

bool TargetDecoyCandidatePairManager::isInit() {
    return !m_targetDecoyCandidatePairs.empty();
}

int TargetDecoyCandidatePairManager::targetsCount() {
    return m_targetDecoyCandidatePairs.size();
}
