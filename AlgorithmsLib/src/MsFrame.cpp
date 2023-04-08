//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "DeisotoperTandem.h"
#include "ErrorUtils.h"
#include "MsScansDenoiseTron.h"

#include <QElapsedTimer>


Err MsFrame::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    m_frame = scanPoints;

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

    QMap<ScanNumber, ScanPoints> deisotopedTandemScans;
    e = DeisotoperTandem::deisotopeTandemScansParallel(
            m_frame,
            m_params.ms2ExtractionWidthPPM,
            &deisotopedTandemScans
    ); ree;

    m_frame = deisotopedTandemScans;

    ERR_RETURN
}

Err MsFrame::smoothFrame() {
    return Error::eFunctionNotImplemented;
}
