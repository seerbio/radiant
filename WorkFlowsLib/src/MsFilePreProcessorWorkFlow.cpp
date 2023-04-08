//
// Created by anichols on 4/7/23.
//

#include "MsFilePreProcessorWorkFlow.h"

#include "DeisotoperTandem.h"
#include "ErrorUtils.h"
#include "MsScansDenoiseTron.h"

#include <QElapsedTimer>


Err MsFilePreProcessorWorkFlow::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}

Err MsFilePreProcessorWorkFlow::preprocessTandemScans(MsReaderPointer *msReaderPointer) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    QMap<ScanNumber, ScanPoints> tandemScanPointsToUpdate;

    const int msLevel = 2;
    e = (*msReaderPointer)->getScanPoints(
            msLevel,
            &tandemScanPointsToUpdate
            ); ree;

    {
        QMap<ScanNumber, ScanPoints> denoisedTandemScans;
        e = denoiseTandemScans(
                tandemScanPointsToUpdate,
                &denoisedTandemScans
        );
        ree;

        tandemScanPointsToUpdate = denoisedTandemScans;
    }

    {
        QMap<ScanNumber, ScanPoints> deisotopedTandemScans;
        e = DeisotoperTandem::deisotopeTandemScansParallel(
                tandemScanPointsToUpdate,
                m_params.ms2ExtractionWidthPPM,
                &deisotopedTandemScans
        ); ree;

        tandemScanPointsToUpdate = deisotopedTandemScans;
    }

    {
        //TODO Put scan mass recalibrator here.
    }

    e = (*msReaderPointer)->updateScanPoints(tandemScanPointsToUpdate); ree

    qDebug() << "File preprocessed in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err MsFilePreProcessorWorkFlow::denoiseTandemScans(
        QMap<ScanNumber, ScanPoints> &tandemScans,
        QMap<ScanNumber, ScanPoints> *denoisedTandemScans
        ) {

    ERR_INIT

    MsScansDenoiseTron denoiser;
    e = denoiser.init(m_params); ree;

    e = denoiser.denoiseScansFrame(
            tandemScans,
            denoisedTandemScans
    ); ree;

    ERR_RETURN
}
