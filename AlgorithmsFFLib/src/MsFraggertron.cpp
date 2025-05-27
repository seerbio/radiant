//
// Created by andrewnichols on 5/23/25.
//

#include "MsFraggertron.h"

#include "TargetDecoyCandidatePair.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN MsFraggertron::Private {
	using rTreeMzDouble = double;
	using rTreeMzHashed = MzHashed;
	using rTreeCoor = bg::model::point<rTreeMzDouble, 1, bg::cs::cartesian>;
	using rTreeSearchBox = bg::model::box<rTreeCoor>;
	using rTreePoint = std::pair<rTreeCoor, rTreeMzHashed>;
	using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

	Private();
	~Private();

	Err init(
		QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPsPntrs
		);

	Err findTargetDecoyCandidates(
		double mz,
		double ppmTolerance,
		QVector<TargetDecoyCandidatePair *> *targetDecoyCandidates
		);

private:
	RTree *m_rTree;
	QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *m_mzHashedVsTDCPPntrs;

};

MsFraggertron::Private::Private()
: m_rTree(Q_NULLPTR)
, m_mzHashedVsTDCPPntrs(Q_NULLPTR)
{}

MsFraggertron::Private::~Private() {
	delete m_rTree;
}

Err MsFraggertron::Private::init(
	QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPsPntrs
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(*mzHashedVsTDCPsPntrs); ree;

	m_mzHashedVsTDCPPntrs = mzHashedVsTDCPsPntrs;

	std::vector<rTreePoint> cloudLoader;
	for (auto it = mzHashedVsTDCPsPntrs->begin(); it != mzHashedVsTDCPsPntrs->end(); it++) {

		const MzHashed mzHashed = it.key();
		const auto mz = MathUtils::unHashDecimal<double>(mzHashed, S_GLOBAL_SETTINGS.HASHING_PRECISION_DDA);
		const QVector<TargetDecoyCandidatePair*> &tdcps = it.value();

		rTreeCoor coor(mz);
		cloudLoader.emplace_back(coor, mzHashed);
	}

	std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
		return l.first.get<0>() < r.first.get<0>();
	});

	delete m_rTree;

	constexpr int maxElements = 16;
	m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsFraggertron initiated";

	ERR_RETURN
}

Err MsFraggertron::Private::findTargetDecoyCandidates(
	double mz,
	double ppmTolerance,
	QVector<TargetDecoyCandidatePair*> *targetDecoyCandidates
	) {

	ERR_INIT

	const double mzTol = MathUtils::calculatePPM(mz, ppmTolerance);
	const double mzMin = mz - mzTol;
	const double mzMax = mz + mzTol;

	const rTreeSearchBox queryBox(
		(rTreeCoor(mzMin)),
		rTreeCoor(mzMax)
	);

	std::vector<rTreePoint> result;
	m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

	QHash<TargetDecoyCandidatePair*, Occurrence> occurrences;
	for (const rTreePoint &p : result) {
		const QVector<TargetDecoyCandidatePair*> prs = m_mzHashedVsTDCPPntrs->value(p.second);
		targetDecoyCandidates->append(prs);
		// for (TargetDecoyCandidatePair *p : prs) {
		// 	occurrences[p]++;
		// }
	}

	ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


MsFraggertron::MsFraggertron()
: d_ptr(QScopedPointer<Private>(new Private))
{}

MsFraggertron::~MsFraggertron() {
}

Err MsFraggertron::init(QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPPntrs) {
	ERR_INIT
	e = ErrorUtils::isNotEmpty(*mzHashedVsTDCPPntrs); ree;
	e = d_ptr->init(mzHashedVsTDCPPntrs); ree;
	ERR_RETURN
}

Err MsFraggertron::findTargetDecoyCandidates(
	double mz,
	double ppmTolerance,
	QVector<TargetDecoyCandidatePair *> *targetDecoyCandidates
	) {

	ERR_INIT

	e = d_ptr->findTargetDecoyCandidates(
			mz,
			ppmTolerance,
			targetDecoyCandidates
			); ree;

	ERR_RETURN
}
