//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "DeisotoperTandem.h"
#include "ErrorUtils.h"
#include "MsScansDenoiseTron.h"

#include <QElapsedTimer>


MsFrame::MsFrame()
: m_collisionEnergy(-1.0)
, m_isoWindowLower(-1.0)
, m_isoWindowUpper(-1.0)
, m_precursorTargetMz(-1.0)
{}

Err MsFrame::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    m_frame = scanPoints;

    e = buildFrameIndexVsScanNumber(); ree;

    ERR_RETURN
}

Err MsFrame::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanPoints,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        double collisionEnergy,
        double precursorTargetMz,
        double isoWindowLower,
        double isoWindowUpper
        ) {

    ERR_INIT

    e = init(
            pythiaParameters,
            scanPoints
            ); ree;

    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;

    const double minVal = 0.0;

    for (const double val : {collisionEnergy, precursorTargetMz, isoWindowLower, isoWindowUpper}) {
        e = ErrorUtils::isAboveThreshold(
                val,
                minVal,
                ErrorUtilsParam::ExcludeThreshold
        ); ree;
    }

    e = buildFrameIndexVsScanNumber(); ree;

    m_collisionEnergy = collisionEnergy;
    m_precursorTargetMz = precursorTargetMz;
    m_isoWindowLower = isoWindowLower;
    m_isoWindowUpper = isoWindowUpper;

    ERR_RETURN
}


Err MsFrame::preprocessMsFrame(
        bool denoise,
        bool deisotope,
        bool smooth
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    if (denoise) {
        e = denoiseFrame(); ree;
    }

    if (deisotope) {
        e = deisotopeFrame(); ree;
    }

    if (smooth) {
        e = smoothFrame(); ree;
    }

    qDebug() << "File preprocessed in" << et.elapsed() << "mSec";

    ERR_RETURN

}

Err MsFrame::denoiseFrame() {

    ERR_INIT

    MsScansDenoiseTron denoiser;
    e = denoiser.init(m_params); ree;

    QMap<ScanNumber, ScanPoints> denoisedFrame;
    e = denoiser.denoiseScansFrame(
            m_frame,
            &denoisedFrame
    ); ree;

    m_frame = denoisedFrame;

    ERR_RETURN
}

Err MsFrame::deisotopeFrame() {

    ERR_INIT

    const bool runParallel = false;

    QMap<ScanNumber, ScanPoints> deisotopedTandemScans;
    e = DeisotoperTandem::deisotopeTandemScansParallel(
            m_frame,
            m_params.ms2ExtractionWidthPPM,
            runParallel,
            &deisotopedTandemScans
    ); ree;

    m_frame = deisotopedTandemScans;

    ERR_RETURN
}

Err MsFrame::smoothFrame() {
    return Error::eFunctionNotImplemented;
}

QPair<double, double> MsFrame::precursorMzTargetStartEnd() const {
    return {m_precursorTargetMz - m_isoWindowLower, m_precursorTargetMz + m_isoWindowUpper};
}

UniqueMsInfoScanKey MsFrame::uniqueMsInfoScanKey() const {
    return m_uniqueMsInfoScanKey;
}

Err MsFrame::buildFrameIndexVsScanNumber() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame); ree;

    for (const ScanNumber &sn : m_frame.keys()) {
        m_frameIndexVsScanNumber.insert(m_frameIndexVsScanNumber.size(), sn);
    }

    ERR_RETURN
}

int MsFrame::scanCount() const {
    return m_frameIndexVsScanNumber.size();
}

QMap<FrameIndex, ScanPoints> MsFrame::frameIndexVsScanPoints() const {

    QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints;

    for (const ScanPoints &sp : m_frame) {
        frameIndexVsScanPoints.insert(frameIndexVsScanPoints.size(), sp);
    }

    return frameIndexVsScanPoints;
}


