//
// Created by anichols on 6/7/24.
//

#include "Ms2IonFraggertronManager.h"

#include "MsCalibratomatic.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "CandidateScores.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN Ms2IonFraggertronManager::Private
{

    struct RTreeMS2IonPoint {
        MS2Ion *ms2Ion = nullptr;
        CandidateScores *candidateScores = nullptr;
    };

    using rTreeScanNumber = float;
    using rTreeIntensity = float;
    using rTreeCoor = bg::model::point<float, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, std::pair<MS2Ion, CandidateScores*>> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(
        const MsCalibratomatic &msCalibratomatic,
        const QVector<CandidateScores*> &candidateScoresPntrs
        );

    Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs);

    Err buildRTreeInput(QVector<rTreePoint> *cloudLoader);

    Err extractMs2Points(
        float mzMin,
        float mzMax,
        float scanTimeMin,
        float scanTimeMax,
        QVector<QPair<MS2Ion, CandidateScores*>> *ms2IonsVsCandidateScoresPntrses
        ) const;

private:

    MsCalibratomatic m_msCalibratomatic;
    QVector<CandidateScores*> m_candidateScores;

    RTree *m_rTree;
    bool m_isInit;

};

Ms2IonFraggertronManager::Private::Private()
: m_rTree(Q_NULLPTR)
, m_isInit(false)
{}

Ms2IonFraggertronManager::Private::~Private() {delete m_rTree;}

Err Ms2IonFraggertronManager::Private::init(
    const MsCalibratomatic &msCalibratomatic,
    const QVector<CandidateScores*> &candidateScoresPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msCalibratomatic.isInitRT()); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

    m_candidateScores = candidateScoresPntrs;
    m_msCalibratomatic = msCalibratomatic;

    QVector<rTreePoint> cloudLoader;
    e = buildRTreeInput(&cloudLoader); ree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    m_isInit = true;

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::buildRTreeInput(QVector<rTreePoint> *cloudLoader) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInitRT()); ree;
    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

    for (CandidateScores *cs : m_candidateScores) {

        float predictedScanTime;
        e = m_msCalibratomatic.predictScanTime(cs->iRt(), &predictedScanTime); ree;

        const QVector<MS2Ion> &ms2Ions = cs->ms2Ions();

        for (const MS2Ion &ms2Ion : ms2Ions) {
            rTreeCoor coor(ms2Ion.mz, predictedScanTime);
            rTreePoint point(coor, {ms2Ion, cs});
            cloudLoader->push_back(point);
        }
    }

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::extractMs2Points(
    float mzMin,
    float mzMax,
    float scanTimeMin,
    float scanTimeMax,
    QVector<QPair<MS2Ion, CandidateScores*>>* ms2IonsVsCandidateScoresPntrses
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;

    const rTreeSearchBox queryBox(
        rTreeCoor(mzMin, scanTimeMin),
        rTreeCoor(mzMax, scanTimeMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    ms2IonsVsCandidateScoresPntrses->reserve(static_cast<int>(result.size()));
    for (const auto &[coordinate, ms2Ion_and_ScorePair] : result) {
        auto &[ms2Ion, candidateScore] = ms2Ion_and_ScorePair;
        ms2IonsVsCandidateScoresPntrses->push_back({ms2Ion, candidateScore});
    }

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) {

    ERR_INIT

    MsCalibratomatic msCalibratomatic;;
    msCalibratomatic.setScanTimeStDev(0.8);

    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

    m_candidateScores = candidateScoresPntrs;
    m_msCalibratomatic = msCalibratomatic;

    QVector<rTreePoint> cloudLoader;
    for (CandidateScores *cs : m_candidateScores) {

        float predictedScanTime = cs->iRt();

        const QVector<MS2Ion> &ms2Ions = cs->ms2Ions();

        for (const MS2Ion &ms2Ion : ms2Ions) {
            rTreeCoor coor(ms2Ion.mz, predictedScanTime);
            rTreePoint point(coor, {ms2Ion, cs});
            cloudLoader.push_back(point);
        }
    }

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    m_isInit = true;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


Ms2IonFraggertronManager::Ms2IonFraggertronManager() : d_ptr(QScopedPointer<Private>(new Private())) {}

Ms2IonFraggertronManager::~Ms2IonFraggertronManager() {}

Err Ms2IonFraggertronManager::init(
        const MsCalibratomatic &msCalibratomatic,
        const QVector<CandidateScores*> &candidateScoresPntrs
        ) const {

    ERR_INIT
    e = d_ptr->init(msCalibratomatic, candidateScoresPntrs); ree;
    ERR_RETURN
}

Err Ms2IonFraggertronManager::extractMs2Points(
    float mzMin,
    float mzMax,
    float scanTimeMin,
    float scanTimeMax,
    QVector<QPair<MS2Ion, CandidateScores*>> *ms2IonsVsCandidateScoresPntrses
    ) {
    ERR_INIT
    e = d_ptr->extractMs2Points(
        mzMin,
        mzMax,
        scanTimeMin,
        scanTimeMax,
        ms2IonsVsCandidateScoresPntrses
        ); ree;
    ERR_RETURN
}

Err Ms2IonFraggertronManager::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) const {
    ERR_INIT
    e = d_ptr->initTesting(candidateScoresPntrs); ree;
    ERR_RETURN
}
