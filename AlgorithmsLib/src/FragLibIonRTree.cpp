//
// Created by anichols on 7/1/23.
//

#include "FragLibIonRTree.h"

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "IsotopicDistributionBuilder.h"

#include <QElapsedTimer>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN FragLibIonRTree::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, int> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QMap<PeptideStringWithMods, MS2IonsSeparated> &fragPreds);

    void insertMs2IonsSeparatedToFragLibIons(
            const QMap<IonIndex, MS2Ion> &ms2IonsSeparated,
            const QString &ionType,
            PeptideId peptideId
    );

    Err addIsotopes();

    Err loadpeptideIdVsFragLibIons();
    Err loadRTree();

    Err buildMzHashedVsFragLibIonFrequencePercentages(
            double ppmTolerance,
            QMap<MzHashed, FrequencyPercent> *mzHashVsFreqPct
            );

    Err addFrequencyPercentagesToFragLibIons(const QMap<MzHashed, FrequencyPercent> &mzHashVsFreqPct);

    Err getFragLibIons(
        double mzMin,
        double mzMax,
        double iRTMin,
        double iRTMax,
        QVector<FragLibIon> *foundFragLibIons
            );

public:
    double m_mzMin;
    double m_mzMax;
    double m_iRTMin;
    double m_iRTMax;

private:

    RTree *m_rTree;
    int m_defaultPrecision;

    QMap<Id, FragLibIon> m_fragLibIons;
    QMap<PeptideId, QVector<FragLibIon>> m_peptideIdVsFragLibIons;
    QMap<PeptideId, PeptideStringWithMods> m_peptideIdVsPeptideStringWithMods;

};


FragLibIonRTree::Private::Private()
: m_defaultPrecision(3)
, m_rTree(Q_NULLPTR)
, m_mzMin(-1.0)
, m_mzMax(-1.0)
, m_iRTMin(-1.0)
, m_iRTMax(-1.0)
{}


FragLibIonRTree::Private::~Private() {
    delete m_rTree;
}

Err FragLibIonRTree::Private::init(const QMap<PeptideStringWithMods, MS2IonsSeparated> &fragPreds) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(fragPreds); ree;

    QHash<PeptideStringWithMods, bool> peptideStringWithModsEntered;

    for (auto it = fragPreds.begin(); it != fragPreds.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const MS2IonsSeparated &ms2IonsSeparated = it.value();

        e = ErrorUtils::isFalse(peptideStringWithModsEntered.contains(peptideStringWithMods)); ree;

        const PeptideId newPeptideId = m_peptideIdVsPeptideStringWithMods.size();

        m_peptideIdVsPeptideStringWithMods.insert(newPeptideId, peptideStringWithMods);
        peptideStringWithModsEntered.insert(peptideStringWithMods, true);

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.bIons,
                S_GLOBAL_SETTINGS.B_IONS,
                newPeptideId
                );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.yIons,
                S_GLOBAL_SETTINGS.Y_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.b2Ions,
                S_GLOBAL_SETTINGS.B2_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.y2Ions,
                S_GLOBAL_SETTINGS.Y2_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.bNH3Ions,
                S_GLOBAL_SETTINGS.B_NH3_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.bH2OIons,
                S_GLOBAL_SETTINGS.B_H2O_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.yNH3Ions,
                S_GLOBAL_SETTINGS.Y_NH3_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.yH2OIons,
                S_GLOBAL_SETTINGS.Y_H2O_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.aIons,
                S_GLOBAL_SETTINGS.A_IONS,
                newPeptideId
        );

        insertMs2IonsSeparatedToFragLibIons(
                ms2IonsSeparated.precursorIons,
                S_GLOBAL_SETTINGS.PRECURSOR_IONS,
                newPeptideId
        );
    }

    e = addIsotopes(); ree

    e = loadpeptideIdVsFragLibIons(); ree
    e = loadRTree(); ree

    qDebug() << "FragLibIon RTree loaded in" << et.elapsed() << "mSec";
    qDebug() << "Peptide Count" << m_peptideIdVsPeptideStringWithMods.size() << fragPreds.size();

    ERR_RETURN
}

void FragLibIonRTree::Private::insertMs2IonsSeparatedToFragLibIons(
        const QMap<IonIndex, MS2Ion> &ms2IonsSeparated,
        const QString &ionType,
        PeptideId peptideId
) {

    for (auto it = ms2IonsSeparated.begin(); it != ms2IonsSeparated.end(); it++) {

        const IonIndex ionIndex = it.key();
        const MS2Ion &ms2Ion = it.value();

        FragLibIon fragLibIon;
        fragLibIon.peptideId = peptideId;
        fragLibIon.mz = ms2Ion.mz;
        fragLibIon.intensity = ms2Ion.intensity;
        fragLibIon.ionIndex = ionIndex;
        fragLibIon.ionType = ionType;
        fragLibIon.iRT = 0.0; //TODO replace w/ predicted iRT when you can predict this.
//        fragLibIon.iMobility = 0.0; //TODO turn this on when using ion mobility and then use 3d rtree
        fragLibIon.charge = ms2Ion.ionLabel.contains("^2") ? 2 : 1;

        m_fragLibIons.insert(m_fragLibIons.size(), fragLibIon);
    }

}

namespace {

    void thresholdIsotopes(QVector<double> *vec) {

        const double isotopePercentThreshold = 0.1;

        const auto terminatorLogic = [isotopePercentThreshold](double pct){
            return pct < isotopePercentThreshold;
        };

        const auto terminator = std::remove_if(vec->begin(), vec->end(), terminatorLogic);

        vec->erase(terminator, vec->end());
    }

    QMap<int, QVector<double>> buildIsotopicDistributionTable() {

        QMap<int, QVector<double>> isoDistTable;

        const int massMax = 10000;
        for (int mass = 0; mass < massMax; mass += 100) {

            QVector<double> isoDist
                = IsotopicDistributionBuilder::getIsotopicDistribution(static_cast<double>(mass));

            thresholdIsotopes(&isoDist);

            isoDistTable.insert(mass, isoDist);
        }

        return isoDistTable;
    }

}//namespace
Err FragLibIonRTree::Private::addIsotopes() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree
    e = ErrorUtils::isNotEmpty(m_peptideIdVsPeptideStringWithMods); ree

    const double fragIonPercentMin = 0.025;

    const QMap<int, QVector<double>> isoDistTable = buildIsotopicDistributionTable();
    e = ErrorUtils::isNotEmpty(isoDistTable); ree

    const int indexMultiplier = 100;

    for (const FragLibIon &fli : m_fragLibIons) {

        const double fragMassApprox = fli.mz * fli.charge;

        const int isoDistTableInd
            = indexMultiplier * static_cast<int>(std::round(fragMassApprox / indexMultiplier));

        const QVector<double> &isoDist = isoDistTable.value(isoDistTableInd);

        for (int isoIndex = 1; isoIndex < isoDist.size(); isoIndex++) {

            const double isotopePercent = isoDist.at(isoIndex);
            const double isoIntensity = fli.intensity * isotopePercent;

            if (isoIntensity < fragIonPercentMin) {
                continue;
            }

            FragLibIon fragLibIonIso = fli;
            fragLibIonIso.isIsotope = true;
            fragLibIonIso.intensity = isoIntensity;
            fragLibIonIso.mz += (isoIndex * S_GLOBAL_SETTINGS.ISO_DIFF) / fli.charge;

            m_fragLibIons.insert(m_fragLibIons.size(), fragLibIonIso);
        }

    }

    ERR_RETURN
}

Err FragLibIonRTree::Private::loadpeptideIdVsFragLibIons() {
    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree

    m_peptideIdVsFragLibIons.clear();

    for (const FragLibIon &fli : m_fragLibIons) {
        m_peptideIdVsFragLibIons[fli.peptideId].push_back(fli);
    }

    ERR_RETURN
}

Err FragLibIonRTree::Private::loadRTree() {

    ERR_INIT

    delete m_rTree;
    m_rTree = nullptr;

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree
    e = ErrorUtils::isNotEmpty(m_peptideIdVsPeptideStringWithMods); ree

    std::vector<rTreePoint> cloudLoader;

    for (auto it = m_fragLibIons.begin(); it != m_fragLibIons.end(); it++) {

        const Id id = it.key();
        const FragLibIon &fragLibIon = it.value();

        rTreeCoor coor(fragLibIon.mz, fragLibIon.iRT);
        cloudLoader.emplace_back(coor, id);
    }

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    m_mzMin = m_rTree->bounds().min_corner().get<0>();
    m_iRTMin = m_rTree->bounds().min_corner().get<1>();
    m_mzMax = m_rTree->bounds().max_corner().get<0>();
    m_iRTMax = m_rTree->bounds().max_corner().get<1>();

    qDebug() << "mz range" << m_mzMin << "to" << m_mzMax;
    qDebug() << "iRT range" << m_iRTMin << "to" << m_iRTMax;
    qDebug() << "FragLibIon RTree size:" << m_rTree->size();

    ERR_RETURN
}

Err FragLibIonRTree::Private::buildMzHashedVsFragLibIonFrequencePercentages(
        double ppmTolerance,
        QMap<MzHashed, FrequencyPercent> *mzHashVsFreqPct
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree;

    const auto fragLibIonsSize = static_cast<double>(m_fragLibIons.size());

    for (const FragLibIon &fli : m_fragLibIons) {

        const int mzHashed = MathUtils::hashDecimal(fli.mz, m_defaultPrecision);

        if (mzHashVsFreqPct->contains(mzHashed)) {
            continue;
        }

        const double mzTol = MathUtils::calculatePPM(fli.mz, ppmTolerance);
        const double mzMin = fli.mz - mzTol;
        const double mzMax = fli.mz + mzTol;

        QVector<FragLibIon> foundFragLibIons;

        e = getFragLibIons(
                mzMin,
                mzMax,
                m_iRTMin,
                m_iRTMax,
                &foundFragLibIons
                ); ree

        const double fragIonFrequencyPct = foundFragLibIons.size() / fragLibIonsSize;
        mzHashVsFreqPct->insert(mzHashed, fragIonFrequencyPct);
    }

    qDebug() << "Fragment frequencies built in:" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err FragLibIonRTree::Private::getFragLibIons(
        double mzMin,
        double mzMax,
        double iRTMin,
        double iRTMax,
        QVector<FragLibIon> *foundFragLibIons
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    const rTreeSearchBox queryBox(
            rTreeCoor(mzMin, iRTMin),
            rTreeCoor(mzMax, iRTMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const rTreePoint &rtp : result) {
        const Id fragLibIonId = rtp.second;
        foundFragLibIons->push_back(m_fragLibIons.value(fragLibIonId));
    }

    ERR_RETURN
}

Err FragLibIonRTree::Private::addFrequencyPercentagesToFragLibIons(const QMap<MzHashed, FrequencyPercent> &mzHashVsFreqPct) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzHashVsFreqPct); ree;
    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree;

    const QList<FrequencyPercent> &freqPcts = mzHashVsFreqPct.values();

    const double frequencyPercentMinLowerer = 10.0;
    const double frequencyPercentMin = *std::min(freqPcts.begin(), freqPcts.end()) / frequencyPercentMinLowerer;

    for (auto it = m_fragLibIons.begin(); it != m_fragLibIons.end(); it++) {

        const FrameIndex frameIndex = it.key();
        FragLibIon newFli = it.value();

        const int mzHashed = MathUtils::hashDecimal(newFli.mz, m_defaultPrecision);
        const double mzFrequencyPercent = mzHashVsFreqPct.value(mzHashed);

        newFli.frequencyPercent = std::max(mzFrequencyPercent, frequencyPercentMin);

        m_fragLibIons[frameIndex] = newFli;
    }

    e = loadpeptideIdVsFragLibIons(); ree
    e = loadRTree(); ree

    ERR_RETURN
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


FragLibIonRTree::FragLibIonRTree() : d_ptr(new Private()) {}

FragLibIonRTree::~FragLibIonRTree() {}

Err FragLibIonRTree::init(const QMap<PeptideStringWithMods, MS2IonsSeparated> &fragPreds) {
    ERR_INIT
    d_ptr->init(fragPreds); ree;
    ERR_RETURN
}

Err FragLibIonRTree::buildMzHashedVsFragLibIonFrequencePercentages(
        double ppmTolerance,
        QMap<MzHashed, FrequencyPercent> *mzHashVsFreqPct
        ) {
    ERR_INIT
    e = d_ptr->buildMzHashedVsFragLibIonFrequencePercentages(
            ppmTolerance,
            mzHashVsFreqPct
            ); ree
    ERR_RETURN
}

Err FragLibIonRTree::getFragLibIons(
        double mzMin,
        double mzMax,
        double iRTMin,
        double iRTMax,
        QVector<FragLibIon> *foundFragLibIons
) {
    ERR_INIT

    e = d_ptr->getFragLibIons(
            mzMin,
            mzMax,
            iRTMin,
            iRTMax,
            foundFragLibIons
            ); ree

    ERR_RETURN
}

Err FragLibIonRTree::getFragLibIons(
        double mzMin,
        double mzMax,
        QVector<FragLibIon> *foundFragLibIons
        ) {

    ERR_INIT

    e = d_ptr->getFragLibIons(
            mzMin,
            mzMax,
            d_ptr->m_iRTMin,
            d_ptr->m_iRTMax,
            foundFragLibIons
    ); ree

    ERR_RETURN
}

Err FragLibIonRTree::addFrequencyPercentagesToFragLibIons(const QMap<MzHashed, FrequencyPercent> &mzHashVsFreqPct) {
    ERR_INIT
    e = d_ptr->addFrequencyPercentagesToFragLibIons(mzHashVsFreqPct); ree
    ERR_RETURN
}
