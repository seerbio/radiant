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
    using rTreePoint = std::pair<rTreeCoor, MS2Frag*> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(
		const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions
    	);

    Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs);

    Err buildRTreeInput(
		const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions,
    	QVector<rTreePoint> *cloudLoader
    	);

    Err extractMs2Points(
        float mzMin,
        float mzMax,
        float iRTMin,
		float iRTMax,
		QVector<MS2Frag*> *tdPeptideFrags
        ) const;

private:

	QVector<MS2Frag> m_ms2Frags;
    RTree *m_rTree;
    bool m_isInit;

};

Ms2IonFraggertronManager::Private::Private()
: m_rTree(Q_NULLPTR)
, m_isInit(false)
{}

Ms2IonFraggertronManager::Private::~Private() {delete m_rTree;}

Err Ms2IonFraggertronManager::Private::init(
	const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2Ions); ree;

    QVector<rTreePoint> cloudLoader;
    e = buildRTreeInput(
    	ms2Ions,
    	&cloudLoader
    	); ree;

	std::sort(
		cloudLoader.begin(),
		cloudLoader.end(),
		[](const rTreePoint &l, const rTreePoint &r) {return l.second->mzVal < r.second->mzVal;}
		);

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    m_isInit = true;

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::buildRTreeInput(
	const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions,
	QVector<rTreePoint> *cloudLoader
	) {

    ERR_INIT

	for (const std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> &tpl : ms2Ions) {
		TargetDecoyCandidatePair *tdcpPntr = std::get<0>(tpl);
		const MS2IonsTarget &ms2Target = std::get<1>(tpl);
		const MS2IonsDecoy &ms2Decoy = std::get<2>(tpl);

		constexpr int decoyDoubler = 2;
		constexpr int top8Ions = 8;
		m_ms2Frags.reserve(ms2Ions.size() * decoyDoubler * top8Ions);

		for (const MS2Ion &ms2Ion : ms2Target.mid(0, std::min(top8Ions, ms2Target.size()))) {

			rTreeCoor coor(ms2Ion.mz, tdcpPntr->iRt());

			MS2Frag ms2Frag;
			ms2Frag.tdcpId = tdcpPntr->id;
			ms2Frag.mzVal = ms2Ion.mz;
			ms2Frag.intensityVal = ms2Ion.intensity;

			m_ms2Frags.push_back(ms2Frag);

			rTreePoint point(coor, &m_ms2Frags.back());
			cloudLoader->push_back(point);
		}

		for (const MS2Ion &ms2Ion : ms2Decoy.mid(0, std::min(top8Ions, ms2Decoy.size()))) {

			rTreeCoor coor(ms2Ion.mz, tdcpPntr->iRt());

			MS2Frag ms2Frag;
			ms2Frag.tdcpId = tdcpPntr->id + 1;
			ms2Frag.mzVal = ms2Ion.mz;
			ms2Frag.intensityVal = ms2Ion.intensity;

			m_ms2Frags.push_back(ms2Frag);

			rTreePoint point(coor, &m_ms2Frags.back());
			cloudLoader->push_back(point);
		}

	}

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::extractMs2Points(
    float mzMin,
    float mzMax,
    float iRTMin,
	float iRTMax,
	QVector<MS2Frag*> *tdPeptideFrags
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    
    const rTreeSearchBox queryBox(
        rTreeCoor(mzMin, iRTMin),
        rTreeCoor(mzMax, iRTMax)
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

Ms2IonFraggertronManager::~Ms2IonFraggertronManager() {}

Err Ms2IonFraggertronManager::init(
	const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions
	) const {

	ERR_INIT

	e = d_ptr->init(ms2Ions); ree;

	ERR_RETURN
}

Err Ms2IonFraggertronManager::extractMs2Points(
    float mzMin,
    float mzMax,
    float iRTMin,
	float iRTMax,
	QVector<MS2Frag*> *tdPeptideFrags
    ) {
    ERR_INIT
    e = d_ptr->extractMs2Points(
        mzMin,
        mzMax,
        iRTMin,
        iRTMax,
        tdPeptideFrags
        ); ree;
    ERR_RETURN
}

Err Ms2IonFraggertronManager::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) const {
    ERR_INIT
    e = d_ptr->initTesting(candidateScoresPntrs); ree;
    ERR_RETURN
}
