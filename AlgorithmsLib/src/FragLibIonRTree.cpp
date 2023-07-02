//
// Created by anichols on 7/1/23.
//

#include "FragLibIonRTree.h"


#include "EigenSparseUtils.h"
#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "MsUtils.h"

#include <QElapsedTimer>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN FragLibIonRTree::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, double> ;
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

    Err loadRTree();

private:

    RTree *m_rTree;
    int m_defaultPrecision;

    QMap<Id, FragLibIon> m_fragLibIons;
    QMap<PeptideId, PeptideStringWithMods> m_peptideIdVsPeptideStringWithMods;

};


FragLibIonRTree::Private::Private()
: m_defaultPrecision(4)
, m_rTree(Q_NULLPTR) {}


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

    e = loadRTree(); ree

    qDebug() << "FragLibIon RTree loaded in" << et.elapsed() << "mSec";

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

Err FragLibIonRTree::Private::addIsotopes() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree
    e = ErrorUtils::isNotEmpty(m_peptideIdVsPeptideStringWithMods); ree

    ERR_RETURN
}

Err FragLibIonRTree::Private::loadRTree() {

    ERR_INIT

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
