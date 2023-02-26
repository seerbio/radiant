//
// Created by anichols on 2/23/23.
//

#include "MsScansDenoiseTron.h"

#include "ErrorUtils.h"


Err MsScansDenoiseTron::init(const PythiaParameters &pythiaParameter) {

    ERR_INIT

    m_pythiaParameters = pythiaParameter;
    m_ffpParams.tolerancePPM = pythiaParameter.ms2ExtractionWidthPPM;
    m_ffpParams.skipScanCount = 2;
    m_ffpParams.minScanCount = 3;

    e = ErrorUtils::isTrue(m_ffpParams.isValid()); ree;
    e = m_featureFinderHillBuilder.init(m_ffpParams);

    ERR_RETURN
}

Err MsScansDenoiseTron::denoiseScansFrame(
        const QMap<ScanNumber, ScanPoints> &scansFrame,
        QMap<ScanNumber, ScanPoints> *scansFrameDenoised
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scansFrame); ree;

    QVector<FeatureFinderHill> featureFinderHills;
    e = m_featureFinderHillBuilder.buildHills(
            scansFrame,
            &featureFinderHills
            ); ree;

    for (const FeatureFinderHill &ffh : featureFinderHills) {

        const QVector<int> &scanNumbers = ffh.scanNumbers();
        const QVector<double> &mzVals = ffh.mzVals();
        const QVector<double> &intensities = ffh.intensities();

        e = ErrorUtils::isNotEmpty(scanNumbers); ree;
        e = ErrorUtils::isEqual(scanNumbers.size(), mzVals.size()); ree;
        e = ErrorUtils::isEqual(scanNumbers.size(), intensities.size()); ree;

        for (int i = 0; i < scanNumbers.size(); i++) {

            const double mz = mzVals.at(i);
            const double intensity = intensities.at(i);
            const int scanNumber = scanNumbers.at(i);

            (*scansFrameDenoised)[scanNumber].push_back({mz, intensity});
        }
    }

    ERR_RETURN
}
