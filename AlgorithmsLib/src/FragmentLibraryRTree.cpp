//
// Created by Drucifer on 7/9/2022.
//

#include "FragmentLibraryRTree.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class Q_DECL_HIDDEN FragmentLibraryRTree::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, int> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<FragLibIon> &fragLibIons,
                 const QPair<int, int> &minMaxCharge,
                 double ppmTolerance,
                 double precursorExtractionWindowThomsons);

    QHash<PeptideId, MZION> getPeptidesTableIds(double mz,
                                           double targetMass,
                                           const QPair<double, double> &targetWindow);

    int size();

private:

    RTree *m_rTree;
    QPair<int, int> m_minMaxCharge;
    double m_ppm;
    double m_precursorExtractionWindowThomsons;
    double m_mzMin;
    double m_mzMax;

};


FragmentLibraryRTree::Private::Private()
: m_ppm(10.0)
, m_precursorExtractionWindowThomsons(1.0)
, m_minMaxCharge({2, 4})
, m_mzMin(-1.0)
, m_mzMax(-1.0)
, m_rTree(nullptr)
{}

FragmentLibraryRTree::Private::~Private() {
    delete m_rTree;
}

QPair<double, double> getMzMinMax(const QVector<FragLibIon> &fragLibIons) {

    const auto logic = [](const FragLibIon &l, const FragLibIon &r ){
        return l.mzFrag < r.mzFrag;
    };

    const auto mzMinMax = std::minmax_element(fragLibIons.begin(), fragLibIons.end(), logic);
    return {mzMinMax.first->mzFrag, mzMinMax.second->mzFrag};
}


Err FragmentLibraryRTree::Private::init(
        const QVector<FragLibIon> &fragLibIons,
        const QPair<int, int> &minMaxCharge,
        double ppmTolerance,
        double precursorExtractionWindowThomsons
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibIons); ree;

    m_minMaxCharge = minMaxCharge;
    m_ppm = ppmTolerance;
    m_precursorExtractionWindowThomsons = precursorExtractionWindowThomsons;

    std::vector<rTreePoint> cloudLouder;
    const auto loadLogic = [&](const FragLibIon &fli) {
        rTreeCoor coor(fli.mzFrag, fli.peptideMass);
        return std::make_pair(coor, fli.peptideId);
    };

    std::transform(fragLibIons.begin(), fragLibIons.end(), back_inserter(cloudLouder), loadLogic);

    const int maxElements = 16;
    m_rTree = new RTree(cloudLouder, bgi::dynamic_quadratic(maxElements));

    const QPair<double, double> mzMinMax = getMzMinMax(fragLibIons);

    m_mzMin = mzMinMax.first;
    m_mzMax = mzMinMax.second;

    ERR_RETURN
}


namespace {

    QVector<double> getAcceptableFilterMasses(
            double precursorTargetMz,
            int chargeMin,
            int chargeMax
            ) {

        QVector<double> filterMasses;

        for (int charge = chargeMin; charge <= chargeMax; ++charge) {

            const double peptideMassChargeAtCharge = BiophysicalCalcs::calculateMassFromThomson(
                    precursorTargetMz,
                    charge,
                    0
                    );

            filterMasses.push_back(peptideMassChargeAtCharge);
        }

        return filterMasses;

    }

}//namespace
QHash<PeptideId, MZION> FragmentLibraryRTree::Private::getPeptidesTableIds(
        double mz,
        double targetMass,
        const QPair<double, double> &targetWindow
) {

    const double ppm = MathUtils::calculatePPM(mz, m_ppm);
    const double mzStart = mz - ppm;
    const double mzEnd = mz + ppm;

    QHash<PeptideId, MZION> peptidesTableIdsVsCoor;

    if (mz < m_mzMin || mz > m_mzMax) {
//        std::cout << "LOG IT MZ OUT OF RANGE " << m_mzMin << " " << mz << " " << m_mzMax << std::endl;
        return peptidesTableIdsVsCoor;
    }

    const QVector<double> filterMasses = getAcceptableFilterMasses(
            targetMass,
            m_minMaxCharge.first,
            m_minMaxCharge.second
    );

    for (double mass : filterMasses) {

        const double startPepMass = mass - targetWindow.first - m_precursorExtractionWindowThomsons;
        const double endPepMass = mass + targetWindow.second + m_precursorExtractionWindowThomsons;

        const rTreeSearchBox queryBox(rTreeCoor(mzStart, startPepMass), rTreeCoor(mzEnd, endPepMass));
        std::vector<rTreePoint> result;
        m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

        for (const rTreePoint &rtp : result) {
            peptidesTableIdsVsCoor.insert(rtp.second, rtp.first.get<0>());
        }

    }

    return peptidesTableIdsVsCoor;
};


int FragmentLibraryRTree::Private::size() {
    return static_cast<int>(m_rTree->size());
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


FragmentLibraryRTree::FragmentLibraryRTree()
: d_ptr(new Private())
, m_key(-1)
{}


FragmentLibraryRTree::~FragmentLibraryRTree() {}


Err FragmentLibraryRTree::init(
        const QVector<FragLibIon> &fragLibIons,
        const QPair<int, int> &minMaxCharge,
        double ppmTolerance,
        double precursorExtractionWindowThomsons
) {

    return d_ptr->init(
            fragLibIons,
            minMaxCharge,
            ppmTolerance,
            precursorExtractionWindowThomsons
            );
}

QHash<PeptideId, MZION> FragmentLibraryRTree::getPeptidesTableIds(
        double mz,
        double targetMass,
        const QPair<double, double> &targetWindow
) {
    return d_ptr->getPeptidesTableIds(mz, targetMass, targetWindow);
}

int FragmentLibraryRTree::size() {
    return d_ptr->size();
}

void FragmentLibraryRTree::setKey(int key) {
    m_key = key;
}

int FragmentLibraryRTree::getKey() {
    return m_key;
}
