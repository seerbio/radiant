//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FrameExtractsReader.h"
#include "MsFrameScoretronProcessormatic.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>


namespace {

    Err setFeatureFinderParams(
            const PythiaParameters &pythiaParameters,
            FeatureFinderParameters *featureFinderParameters
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        featureFinderParameters->sigma = pythiaParameters.sigma;
        featureFinderParameters->filterLength = pythiaParameters.filterLength;
        featureFinderParameters->smoothCount = pythiaParameters.smoothCount;
        featureFinderParameters->signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        featureFinderParameters->tolerancePPM = pythiaParameters.ms2ExtractionWidthPPM;
        featureFinderParameters->skipScanCount = pythiaParameters.skipScanCount;
        featureFinderParameters->minScanCount = pythiaParameters.minScanCount;

        featureFinderParameters->printParams();

        e = ErrorUtils::isTrue(featureFinderParameters->isValid()); ree

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    params.print();

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isAboveThreshold(
            mzTargetStartStop.second,
            mzTargetStartStop.first,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;

    m_params = params;
    m_msDataFilePath = msDataFilePath;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_mzTargetStartStop = mzTargetStartStop;

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibFilePath,
            m_params.chargeStateMin,
            m_params.chargeStateMax,
            m_mzTargetStartStop.first,
            m_mzTargetStartStop.second,
            &m_fragPreds,
            &m_fragPredsIsDecoy
            ); ree

    e = ErrorUtils::isNotEmpty(m_fragPreds); ree

    e = initMsFrame(
        msDataFilePath,
        uniqueMsInfoScanKey,
        mzTargetStartStop
        ); ree

    qDebug() << "TargetKey" << uniqueMsInfoScanKey;
    qDebug() << "Scan Count" << m_msFrame.scanCount();
    qDebug() << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}



Err MsFrameScoretron::initMsFrame(
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    e = MsFrame::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop,
            &m_msFrame
    ); ree;

    e = m_msFrame.filterFrameByMz(
            m_params.mzMinDataStructure,
            m_params.mzMaxDataStructure
    ); ree;

//    e = m_msFrame.smoothFrame(
//            m_params.filterLength,
//            m_params.sigma,
//            m_params.smoothCount,
//            m_params.mzMaxDataStructure
//    ); ree

    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
    );ree;




    ERR_RETURN
}

namespace {

    Err buildMsFrameApexTurboXIC(
            const PythiaParameters &pythiaParameters,
            const QMap<FrameIndex, ScanPoints> &frameIndexVsScanPoints,
            TurboXIC *turboXic
            ) {

        ERR_INIT

        FeatureFinderParameters featureFinderParameters;
        e = setFeatureFinderParams(pythiaParameters, &featureFinderParameters); ree

        FeatureFinderHillBuilder featureFinderHillBuilder;
        e = featureFinderHillBuilder.init(featureFinderParameters); ree

        e = featureFinderHillBuilder.buildHills(frameIndexVsScanPoints); ree

        QVector<FeatureFinderHill> featureFinderHills;
        e = featureFinderHillBuilder.featureFinderHills(&featureFinderHills); ree

        //TODO consider making these settable in parameters
        const int smoothCount = 2;
        const int filterLength = 5;
        const int order = 1;
        const int derivative = 0;
        const int rate = 1;

        QMap<FrameIndex, ScanPoints> frameApexes;

        for (const FeatureFinderHill &ffh : featureFinderHills) {

            Eigen::VectorX<double> intensitiesVec = EigenUtils::convertQVectorToEigenVector(ffh.intensities());

            for (int i = 0; i < smoothCount; i++) {
                e = EigenKernelUtils::savitskyGolaySmooth(
                        filterLength,
                        order,
                        derivative,
                        rate,
                        &intensitiesVec
                        ); ree
            }

            const QVector<FrameIndex> &frameIndexes = ffh.scanNumbers();
            const double mzMean = ffh.mzMean();
            const QVector<double> &originalIntensities = ffh.intensities();

            const QMap<int, double> apexes = EigenUtils::apexes(intensitiesVec, 1);

            for (auto it = apexes.begin(); it != apexes.end(); it++) {

                const int apexIndex = it.key();
                const double _smoothedIntensityApexValue = it.value();
                const int frameIndex = frameIndexes.at(apexIndex);
                const double originalIntensity = originalIntensities.at(apexIndex);

                frameApexes[frameIndex].push_back({mzMean, originalIntensity});
            }

        }

        turboXic->init(frameApexes); ree

        ERR_RETURN

    }

    Err buildApexSpectra(
            const TurboXIC &frameApexesTurboXIC,
            QMap<FrameIndex, ScanPoints> *apexScanPoints
            ) {

        ERR_INIT

        double frameIndexMin;
        double frameIndexMax;
        double mzMin;
        double mzMax;

        e = frameApexesTurboXIC.getRTreeLimits(
                &frameIndexMin,
                &frameIndexMax,
                &mzMin,
                &mzMax
                ); ree

        for (
            auto frameIndex = static_cast<int>(frameIndexMin);
            frameIndex <= static_cast<int>(frameIndexMax);
            ++frameIndex
            ) {

            const int frameIndexLo = std::max(0, frameIndex - 1);
            const int frameIndexHi = std::min(static_cast<int>(frameIndexMax), frameIndex + 1);

            const ScanPoints frameIndexScanPoints = frameApexesTurboXIC.extractSpectrum(
                    mzMin,
                    mzMax,
                    frameIndexLo,
                    frameIndexHi
                    );

            apexScanPoints->insert(frameIndex, frameIndexScanPoints);
        }

        ERR_RETURN
    }


}//namespace
Err MsFrameScoretron::scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates) {

    ERR_INIT

    TurboXIC frameApexesTurboXIC;

    e = buildMsFrameApexTurboXIC(
            m_params,
            m_msFrame.frameIndexVsScanPoints(),
            &frameApexesTurboXIC
            ); ree

    QMap<FrameIndex, ScanPoints> apexScanPoints;
    e = buildApexSpectra(
            frameApexesTurboXIC,
            &apexScanPoints
            ); ree

//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    for(auto it = apexScanPoints.begin(); it != apexScanPoints.end(); it++) {
        const ScanNumber sn = m_msFrame.scanNumberFromFrameIndex(it.key());
        scanNumberVsScanPoints.insert(sn, it.value());
    }

    e = MsFrame::writeFrameScans(scanNumberVsScanPoints, "frame.prq"); ree
#endif

    ERR_RETURN
}
