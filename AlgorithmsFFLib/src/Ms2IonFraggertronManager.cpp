//
// Created by anichols on 6/7/24.
//

#include "Ms2IonFraggertronManager.h"

#include "MsCalibratomatic.h"
#include "TargetDecoyCandidatePair.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>


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
    using rTreePoint = std::pair<rTreeCoor, MS2IonLibrary*> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<MS2IonLibrary*> &ms2IonsLibrary);

    Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs);

    Err buildRTreeInput(QVector<rTreePoint> *cloudLoader);

    Err extractMs2Points(
	    float mzPrecursorMin,
		float mzPrecursorMax,
		float mzMin,
		float mzMax,
		QVector<MS2IonLibrary*> *tdPeptideFrags
        ) const;

private:

	QVector<MS2IonLibrary*> m_ms2IonsLibrary;
    RTree *m_rTree;
    bool m_isInit;

};

Ms2IonFraggertronManager::Private::Private()
: m_rTree(Q_NULLPTR)
, m_isInit(false)
{}

Ms2IonFraggertronManager::Private::~Private() {delete m_rTree;}

Err Ms2IonFraggertronManager::Private::init(
	const QVector<MS2IonLibrary*> &ms2IonsLibrary
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsLibrary); ree;
	m_ms2IonsLibrary = ms2IonsLibrary;

    QVector<rTreePoint> cloudLoader;
    e = buildRTreeInput(&cloudLoader); ree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    m_isInit = true;

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::buildRTreeInput(QVector<rTreePoint> *cloudLoader) {

    ERR_INIT

	e = ErrorUtils::isNotEmpty(m_ms2IonsLibrary); ree;

	for (MS2IonLibrary *msil : m_ms2IonsLibrary) {
		rTreeCoor coor(msil->targeDecoyCandidatePairPntr->mz(msil->isDecoy), msil->ms2IonPntr->mz);
		rTreePoint point(coor, msil);
		cloudLoader->push_back(point);
	}

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::extractMs2Points(
	float mzPrecursorMin,
	float mzPrecursorMax,
	float mzMin,
	float mzMax,
	QVector<MS2IonLibrary*> *tdPeptideFrags
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    
    const rTreeSearchBox queryBox(
        rTreeCoor(mzPrecursorMin, mzMin),
        rTreeCoor(mzPrecursorMax, mzMax)
    );
    
    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    tdPeptideFrags->reserve(static_cast<int>(result.size()));
    for (const auto &[fst, snd] : result) {
        tdPeptideFrags->push_back(snd);
    }

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) {

    ERR_INIT

    MsCalibratomatic msCalibratomatic;;
    msCalibratomatic.setScanTimeStDev(0.8);

    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

    // m_candidateScores = candidateScoresPntrs;
    // m_msCalibratomatic = msCalibratomatic;
    //
    // QVector<rTreePoint> cloudLoader;
    // for (CandidateScores *cs : m_candidateScores) {
    //
    //     float predictedScanTime = cs->targetDecoyCandidatePair->iRt();
    //
    //     const QVector<MS2Ion> &ms2Ions = cs->isDecoy
    //                                    ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
    //                                    : cs->targetDecoyCandidatePair->ms2IonsTarget();
    //
    //     for (const MS2Ion &ms2Ion : ms2Ions) {
    //         rTreeCoor coor(ms2Ion.mz, predictedScanTime);
    //         rTreePoint point(coor, {ms2Ion, cs});
    //         cloudLoader.push_back(point);
    //     }
    // }
    //
    // constexpr int maxElements = 16;
    // m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    //
    // m_isInit = true;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


Ms2IonFraggertronManager::Ms2IonFraggertronManager() : d_ptr(QScopedPointer<Private>(new Private())) {}

Ms2IonFraggertronManager::~Ms2IonFraggertronManager() = default;

Err Ms2IonFraggertronManager::init(
	const QVector<MS2IonLibrary*> &ms2IonsLibrary
	) const {

	ERR_INIT
	e = d_ptr->init(ms2IonsLibrary); ree;
	ERR_RETURN
}

Err Ms2IonFraggertronManager::extractMs2Points(
	float mzPrecursorMin,
	float mzPrecursorMax,
    float mzMin,
    float mzMax,
	QVector<MS2IonLibrary*> *tdPeptideFrags
    ) const {
    ERR_INIT
    e = d_ptr->extractMs2Points(
    	mzPrecursorMin,
    	mzPrecursorMax,
        mzMin,
        mzMax,
        tdPeptideFrags
        ); ree;
    ERR_RETURN
}

Err Ms2IonFraggertronManager::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) const {
    ERR_INIT
    e = d_ptr->initTesting(candidateScoresPntrs); ree;
    ERR_RETURN
}
