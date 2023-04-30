//
// Created by anichols on 4/29/23.
//

#include "Ms1FeatureFinder.h"

#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsFrame.h"
#include "MsReaderParquet.h"

Err Ms1FeatureFinder::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}

namespace {

    Err buildMsFrame(
            const QString &msDataFilePath,
            const PythiaParameters &params,
            MsFrame *msFrame,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
            ) {

        ERR_INIT

        MsReaderParquet msReader;
        e = msReader.openFile(
                msDataFilePath,
                "msLevel",
                {1, 1}
        ); ree;

        const int msLevel = 1;
        QMap<ScanNumber, ScanPoints> ms1Scans;
        e = msReader.getScanPoints(
                msLevel,
                &ms1Scans
        ); ree;

        e = msFrame->init(
                params,
                "-1000",
                ms1Scans,
                {-1.0, -1.0}
        ); ree;



        const bool denoise = true;
        const bool deisotope = false;
        const bool smooth = true;
        e = msFrame->preprocessMsFrame(
                denoise,
                deisotope,
                smooth
        ); ree;

        *scanNumberVsScanTime = msReader.getScanNumberVsScanTime();

        e = ErrorUtils::isNotEmpty(*scanNumberVsScanTime); ree;

        ERR_RETURN
    }


}//namespace
Err Ms1FeatureFinder::exec(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsFrame msFrame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    e = buildMsFrame(
            msDataFilePath,
            m_params,
            &msFrame,
            &scanNumberVsScanTime
    ); ree;

    //TODO set these accordingly in PythiaParams file.
    FeatureFinderParameters ffParams;
    ffParams.tolerancePPM = m_params.featureFinderTolerancePPM;
    ffParams.skipScanCount = 0;
    ffParams.minScanCount = 3;
    ffParams.useMeanMz = true;

    FeatureFinderHillBuilder hillBuilder;
    e = hillBuilder.init(ffParams); ree;
    hillBuilder.setRunParallel(true);

    QVector<FeatureFinderHill> featureFinderHills;
    e = hillBuilder.buildHills(
            msFrame.scanNumberVsScanPoints(),
            &featureFinderHills
            ); ree;

    qDebug() << featureFinderHills.size();
    e = hillBuilder.refineHills(&featureFinderHills); ree;
    qDebug() << featureFinderHills.size();

#define WRITE_HILLS
#ifdef WRITE_HILLS
    const QString batmassHillsFilePath = msDataFilePath + ".newRefine" + ".mzrt.csv";
    e = FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
            scanNumberVsScanTime,
            featureFinderHills,
            batmassHillsFilePath
            ); ree;
#endif

    ERR_RETURN
}
