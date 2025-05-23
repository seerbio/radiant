//
// Created by andrewnichols on 5/23/25.
//

#include "MsFraggatron.h"

#include "TargetDecoyCandidatePair.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN MsFraggatron::Private {
	using rTreeMzDouble = double;
	using rTreeMzHashed = MzHashed;
	using rTreeCoor = bg::model::point<rTreeMzDouble, 1, bg::cs::cartesian>;
	using rTreeSearchBox = bg::model::box<rTreeCoor>;
	using rTreePoint = std::pair<rTreeCoor, rTreeMzHashed>;
	using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

	Private();
	~Private();

	Err init(QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPsPntrs);

private:
	RTree *m_rTree;

};

MsFraggatron::Private::Private()
: m_rTree(Q_NULLPTR)
{}

MsFraggatron::Private::~Private() {
	delete m_rTree;
}

Err MsFraggatron::Private::init(QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPsPntrs) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(*mzHashedVsTDCPsPntrs); ree;

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

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsFraggatron initiated";

	ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


MsFraggatron::MsFraggatron()
: m_mzHashedVsTDCPPntrs(Q_NULLPTR)
, d_ptr(QScopedPointer<Private>(new Private))
{}

MsFraggatron::~MsFraggatron() {
}

Err MsFraggatron::init(QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPPntrs) {
	ERR_INIT
	e = ErrorUtils::isNotEmpty(*mzHashedVsTDCPPntrs); ree;
	m_mzHashedVsTDCPPntrs = mzHashedVsTDCPPntrs;
	e = d_ptr->init(mzHashedVsTDCPPntrs); ree;
	ERR_RETURN
}
