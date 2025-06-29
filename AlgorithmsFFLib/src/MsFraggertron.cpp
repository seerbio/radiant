//
// Created by andrewnichols on 5/23/25.
//

#include "MsFraggertron.h"

#include "MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePair.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN MsFraggertron::Private {
	using rTreeMzDouble = double;
	using rTreeScanPointsPntrs = ScanPoints;
	using rTreeCoor = bg::model::point<rTreeMzDouble, 1, bg::cs::cartesian>;
	using rTreeSearchBox = bg::model::box<rTreeCoor>;
	using rTreePoint = std::pair<rTreeSearchBox, rTreeScanPointsPntrs*>;
	using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

	Private();
	~Private();

	Err init(MsReaderPointerAcc *msReader);

	Err findScanPointsPntrs(
		double mz,
		double precursorExtractionWindowThomsons,
		QVector<ScanPoints*> *scanPointsPntrs
		) const;

private:
	RTree *m_rTree;
	QMap<ScanNumber, MsScanInfo> m_scanInfos;
	QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *m_mzHashedVsTDCPPntrs;

};

MsFraggertron::Private::Private()
: m_rTree(Q_NULLPTR)
, m_mzHashedVsTDCPPntrs(Q_NULLPTR)
{}

MsFraggertron::Private::~Private() {
	delete m_rTree;
}

Err MsFraggertron::Private::init(MsReaderPointerAcc *msReader) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReader->isInit()); ree;

	constexpr int msLevel = 2;
	m_scanInfos = msReader->ptr->getMsScanInfos(msLevel);
	QMap<ScanNumber, ScanPoints*> scanPointsPntrs;
	e = msReader->ptr->getScanPoints(msLevel, &scanPointsPntrs);

	std::vector<rTreePoint> cloudLoader;
	for (auto it = m_scanInfos.begin(); it != m_scanInfos.end(); it++) {

		const ScanNumber scanNumber = it.key();
		const MsScanInfo &msi = it.value();
		ScanPoints* sp = scanPointsPntrs.value(scanNumber);


		const double mzMin = msi.precursorTargetMz - msi.isoWindowLower;
		const double mzMax = msi.precursorTargetMz + msi.isoWindowLower;

		const rTreeSearchBox windowBox(
			(rTreeCoor(mzMin)),
			(rTreeCoor(mzMax))
		);

		cloudLoader.emplace_back(windowBox, sp);
	}

	delete m_rTree;

	constexpr int maxElements = 16;
	m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsFraggertron initiated";

	ERR_RETURN
}

Err MsFraggertron::Private::findScanPointsPntrs(
	double mz,
	double precursorExtractionWindowThomsons,
	QVector<ScanPoints*> *scanPointsPntrs
	) const {

	ERR_INIT

	const double mzMin = mz - precursorExtractionWindowThomsons;
	const double mzMax = mz + precursorExtractionWindowThomsons;

	const rTreeSearchBox queryBox(
		(rTreeCoor(mzMin)),
		(rTreeCoor(mzMax))
	);

	std::vector<rTreePoint> result;
	m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

	for (const rTreePoint &p : result) {
		scanPointsPntrs->push_back(p.second);
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

Err MsFraggertron::init(MsReaderPointerAcc *msReader) {
	ERR_INIT
	e = ErrorUtils::isTrue(msReader->isInit()); ree;
	e = d_ptr->init(msReader); ree;
	ERR_RETURN
}

Err MsFraggertron::findScanPointsPntrs(
	double mz,
	double precursorExtractionWindowThomsons,
	QVector<ScanPoints*> *scanPointsPntrs
	) {

	ERR_INIT

	e = d_ptr->findScanPointsPntrs(
			mz,
			precursorExtractionWindowThomsons,
			scanPointsPntrs
			); ree;

	ERR_RETURN
}
