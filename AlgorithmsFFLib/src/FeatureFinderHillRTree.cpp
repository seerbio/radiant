//
// Created by anichols on 12/8/23.
//

#include "FeatureFinderHillRTree.h"

#include "FeatureFinderHill.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN FeatureFinderHillRTree::Private
{
    using RTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using RTreeBox = bg::model::box<RTreeCoor>;
    using RTreePoint = std::pair<RTreeBox, FeatureFinderHill*> ;
    using RTree = bgi::rtree<RTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<FeatureFinderHill*> &featureFinderHills);

    Err getRTreeLimits(
            double *mzMin,
            double *mzMax,
            FrameIndex *frameIndexMin,
            FrameIndex *FrameIndexMax
    );

    bool isInit();

    Err getHills(
            double mzMin,
            double mzMax,
            FrameIndex frameIndexMin,
            FrameIndex FrameIndexMax,
            QVector<FeatureFinderHill*> *featureFinderHills
    );

private:

    RTree *m_rTree;

};


FeatureFinderHillRTree::Private::Private() : m_rTree(Q_NULLPTR) {}

FeatureFinderHillRTree::Private::~Private() {
    delete m_rTree;
}

Err FeatureFinderHillRTree::Private::init(const QVector<FeatureFinderHill*> &featureFinderHills) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

    std::vector<RTreePoint> cloudLoader;
    cloudLoader.reserve(featureFinderHills.size());
    std::transform(
            featureFinderHills.begin(),
            featureFinderHills.end(),
            std::back_inserter(cloudLoader),
            [](FeatureFinderHill* ffh){

                const QPair<double, double> mzMinMax = ffh->mzMinMax();
                const QPair<ScanNumber, ScanNumber> intensityMinMax = ffh->scanNumberIndexMinMax();

                const RTreeBox box(
                        RTreeCoor(mzMinMax.first, static_cast<double>(intensityMinMax.first)),
                        RTreeCoor(mzMinMax.second, static_cast<double>(intensityMinMax.second))
                        );
                const RTreePoint rTreePoint(box, ffh);
                return rTreePoint;
            }
            );

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}

Err FeatureFinderHillRTree::Private::getRTreeLimits(
        double *mzMin,
        double *mzMax,
        FrameIndex *frameIndexMin,
        FrameIndex *FrameIndexMax
) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *mzMin = m_rTree->bounds().min_corner().get<0>();
    *mzMax = m_rTree->bounds().max_corner().get<0>();
    *frameIndexMin = static_cast<FrameIndex>(m_rTree->bounds().min_corner().get<1>());
    *FrameIndexMax = static_cast<FrameIndex>(m_rTree->bounds().max_corner().get<1>());

    ERR_RETURN
}

bool FeatureFinderHillRTree::Private::isInit() {
    return !m_rTree->empty();
}

Err FeatureFinderHillRTree::Private::getHills(
        double mzMin,
        double mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        QVector<FeatureFinderHill*> *featureFinderHills
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(isInit()); ree;

    const RTreeBox queryBox(
        RTreeCoor(mzMin, static_cast<double>(frameIndexMin)),
        RTreeCoor(mzMax,static_cast<double>(frameIndexMax))
    );

    std::vector<RTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const RTreePoint &rtp : result) {
        featureFinderHills->push_back(rtp.second);
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

FeatureFinderHillRTree::FeatureFinderHillRTree() : d_ptr(new Private()) {}

FeatureFinderHillRTree::~FeatureFinderHillRTree() {}

Err FeatureFinderHillRTree::init(const QVector<FeatureFinderHill*> &featureFinderHills) {
    ERR_INIT
    e = d_ptr->init(featureFinderHills); ree;
    ERR_RETURN
}

Err FeatureFinderHillRTree::getHills(
        double mzMin,
        double mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        QVector<FeatureFinderHill*> *featureFinderHills
        ) {
    ERR_INIT
    e = d_ptr->getHills(
            mzMin,
            mzMax,
            frameIndexMin,
            frameIndexMax,
            featureFinderHills
            ); ree;
    ERR_RETURN
}
