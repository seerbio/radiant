//
// Created by anichols on 6/4/24.
//

#include "SpectrumCentricMzTargetFrameSearch.h"

#include "MsUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using rTreeTargetDecoyPtr = TargetDecoyCandidatePair*;
using rTreeCoor = bg::model::point<float, 1, bg::cs::cartesian>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeCoor, rTreeTargetDecoyPtr>;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

Err SpectrumCentricMzTargetFrameSearch::init(
    const PythiaParameters& pythiaParameters,
    const MsCalibratomatic& msCalibratomatic,
    const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
    const QVector<TargetDecoyCandidatePair*>& targetDecoyCandidatePairs,
    const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(diaTargetFrame); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msCalibratomatic.isInitRT()); ree;

    m_msCalibratomatic = msCalibratomatic;
    m_pythiaParameters = pythiaParameters;
    m_diaTargetFrame = diaTargetFrame;
    m_targetDecoyCandidatePairs = targetDecoyCandidatePairs;
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    ERR_RETURN
}

namespace {

    Err buildRTreeInput(
        const MsCalibratomatic &msCalibratomatic,
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        QVector<rTreePoint> *cloudLoader
        ) {

        ERR_INIT
        e = ErrorUtils::isTrue(msCalibratomatic.isInitRT()); ree;
        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairs); ree;

        for (TargetDecoyCandidatePair* tdcp : targetDecoyCandidatePairs) {
            float predictedScanTime;
            e = msCalibratomatic.predictScanTime(tdcp->iRt(), &predictedScanTime); ree;

            rTreeCoor coor(predictedScanTime);
            rTreePoint point(coor, tdcp);
            cloudLoader->push_back(point);
        }

        ERR_RETURN
    }

    QVector<QPointF> ms2IonsToQPoint(const QVector<MS2Ion> &ms2Ions) {

        QVector<QPointF> qpoints;
        std::transform(
            ms2Ions.begin(),
            ms2Ions.end(),
            std::back_inserter(qpoints),
            [](const MS2Ion &ms2Ion){return QPointF(ms2Ion.mz, ms2Ion.intensity);}
            );

        return qpoints;
    }

}
Err SpectrumCentricMzTargetFrameSearch::assignIdsToScans() {

    ERR_INIT

    QVector<rTreePoint> cloudLoader;
    e = buildRTreeInput(
        m_msCalibratomatic,
        m_targetDecoyCandidatePairs,
        &cloudLoader
        ); ree;

    constexpr int maxElements = 16;
    RTree rTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    const float scanTimeWindow = m_msCalibratomatic.scanTimeStDev(S_GLOBAL_SETTINGS.STDEV_MULTIPLIER);

    for (auto it = m_diaTargetFrame.begin(); it != m_diaTargetFrame.end(); ++it) {

        const ScanNumber scanNumber = it.key();
        const ScanTime scanTime = m_scanNumberVsScanTime.value(scanNumber);
        const ScanPoints* scanPoints = it.value();

        QVector<QPointF> scanPointsDouble;
        std::transform(
            scanPoints->begin(),
            scanPoints->end(),
            std::back_inserter(scanPointsDouble),
            [](const ScanPoint &sp){return QPointF(sp.x(), sp.y());}
            );

        const ScanTime scanTimeMin = scanTime - scanTimeWindow;
        const ScanTime scanTimeMax = scanTime + scanTimeWindow;

        const rTreeSearchBox queryBox(
            (rTreeCoor(scanTimeMin)),
            rTreeCoor(scanTimeMax)
            );

        std::vector<rTreePoint> result;
        rTree.query(bgi::intersects(queryBox), std::back_inserter(result));

        for (const rTreePoint &p : result) {
            const QVector<MS2Ion> &targetIons = p.second->ms2IonsTarget();
            const QVector<QPointF> &targetIonsPoints = ms2IonsToQPoint(targetIons);
            const QVector<MS2Ion> &decoyIons = p.second->ms2IonsDecoy();
            const QVector<QPointF> &decoyIonsPoints = ms2IonsToQPoint(decoyIons);

            const ExtractPoints epTarget = MsUtils::extractPointsFromPoints(
                scanPointsDouble,
                targetIonsPoints,
                m_pythiaParameters.ms2ExtractionWidthPPM
                );
            const ExtractPoints epDecoy = MsUtils::extractPointsFromPoints(
                scanPointsDouble,
                decoyIonsPoints,
                m_pythiaParameters.ms2ExtractionWidthPPM
                );
        }

    }

    ERR_RETURN
}
