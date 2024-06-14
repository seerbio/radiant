//
// Created by anichols on 6/4/24.
//

#include "SpectrumCentricMzTargetFrameSearch.h"

#include "Deconvolvotron.h"
#include "Ms2IonFraggertronManager.h"
#include "MsUtils.h"


Err SpectrumCentricMzTargetFrameSearch::init(
    const PythiaParameters& pythiaParameters,
    const MsCalibratomatic& msCalibratomatic,
    const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
    const QVector<CandidateScores*> &candidateScoresPntrs,
    const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(diaTargetFrame); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msCalibratomatic.isInitRT()); ree;

    m_msCalibratomatic = msCalibratomatic;
    m_pythiaParameters = pythiaParameters;
    m_diaTargetFrame = diaTargetFrame;
    m_candidateScoresPntrs = candidateScoresPntrs;
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    ERR_RETURN
}

namespace {

    void filterFoundCandidates(QHash<CandidateScores*, QVector<MS2Ion>> *candidateVsMs2IonsFound) {

        QHash<CandidateScores*, QVector<MS2Ion>> candidateVsMs2IonsFoundFiltered;
        for (auto it = candidateVsMs2IonsFound->begin(); it != candidateVsMs2IonsFound->end(); ++it) {

            CandidateScores *cs = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value();

            constexpr int minFoundPointsSize = 3;
            if (ms2Ions.size() < minFoundPointsSize) {
                continue;
            }

            // constexpr ushort rankCountMin = 3;
            // const long top3Ms2IonsCount = std::count_if(
            //     ms2Ions.begin(),
            //     ms2Ions.end(),
            //     [rankCountMin](const MS2Ion &ms2Ion){return ms2Ion.rank <= rankCountMin;}
            //     );
            //
            // constexpr int rankCountFoundMin = 2;
            // if (top3Ms2IonsCount < rankCountFoundMin) {
            //     continue;
            // }

            candidateVsMs2IonsFoundFiltered.insert(cs, ms2Ions);
        }

        *candidateVsMs2IonsFound = candidateVsMs2IonsFoundFiltered;
    }

}//namespace
Err SpectrumCentricMzTargetFrameSearch::assignIdsToScans(
    QVector<QPair<CandidateScores*, DeconvolvotronResult>> *candidateScoresPntrVsScoreses
    ) {

    ERR_INIT

    const float scanTimeWindow
        = m_msCalibratomatic.scanTimeStDev(static_cast<float>(S_GLOBAL_SETTINGS.STDEV_MULTIPLIER));

    constexpr int mzGroupingPrecisioin = 1;
    Deconvolvotron deconvolvotron;
    e = deconvolvotron.init(
        mzGroupingPrecisioin,
        m_pythiaParameters.ms2ExtractionWidthPPM
        ); ree;

    Ms2IonFraggertronManager ms2IonFraggertronManager;
    e = ms2IonFraggertronManager.init(
        m_msCalibratomatic,
        m_candidateScoresPntrs
        ); ree;

    FrameIndex frameIndex = 0;
    for (auto it = m_diaTargetFrame.begin(); it != m_diaTargetFrame.end(); ++it) {

        frameIndex++;

        const ScanNumber scanNumber = it.key();
        const ScanTime scanTime = m_scanNumberVsScanTime.value(scanNumber);
        const ScanPoints *scanPoints = it.value();

        const ScanTime scanTimeMin = scanTime - scanTimeWindow;
        const ScanTime scanTimeMax = scanTime + scanTimeWindow;

        QHash<CandidateScores*, QVector<MS2Ion>> candidateVsMs2IonsFound;

        for (const ScanPoint &scanPoint : *scanPoints) {

            const float massTol = MathUtils::calculatePPM(
                scanPoint.x(),
                static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
                );

            const float mzMin = scanPoint.x() - massTol;
            const float mzMax = scanPoint.x() + massTol;

            QVector<QPair<MS2Ion, CandidateScores*>> ms2IonsVsCandidateScoresPntrses;
            e = ms2IonFraggertronManager.extractMs2Points(
                mzMin,
                mzMax,
                scanTimeMin,
                scanTimeMax,
                &ms2IonsVsCandidateScoresPntrses
                ); ree;

            for (const QPair<MS2Ion, CandidateScores*> &pr : ms2IonsVsCandidateScoresPntrses) {
                candidateVsMs2IonsFound[pr.second].push_back(pr.first);
            }

        }

        if (candidateVsMs2IonsFound.isEmpty()) {
            continue;
        }

        filterFoundCandidates(&candidateVsMs2IonsFound);

        if (candidateVsMs2IonsFound.isEmpty()) {
            continue;
        }

        QVector<QPair<CandidateScores*, QVector<QPointF>>> aMatrixPoints;
        aMatrixPoints.reserve(candidateVsMs2IonsFound.size() * 2);
        for (auto it = candidateVsMs2IonsFound.begin(); it != candidateVsMs2IonsFound.end(); ++it) {
            CandidateScores *cs = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value();

            QVector<QPointF> points;
            std::transform(
                ms2Ions.begin(),
                ms2Ions.end(),
                std::back_inserter(points),
                [](const MS2Ion &ms2Ion){return QPointF(ms2Ion.mz, ms2Ion.intensity);}
                );

            aMatrixPoints.push_back({cs, points});
        }

        QVector<QPointF> scanPointsDouble;
        std::transform(
            scanPoints->begin(),
            scanPoints->end(),
            std::back_inserter(scanPointsDouble),
            [](const ScanPoint &scanPoint){return QPointF(scanPoint.x(), scanPoint.y());}
            );

        QVector<QPair<CandidateScores*, DeconvolvotronResult>> candidateScoresPntrVsScores;
        e = deconvolvotron.deconvolve(
            aMatrixPoints,
            scanPointsDouble,
            &candidateScoresPntrVsScores
            ); ree;
        candidateScoresPntrVsScoreses->append(candidateScoresPntrVsScores);
    }

    ERR_RETURN
}
