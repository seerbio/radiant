//
// Created by Drucifer on 12/4/2021.
//

#include "PythiaParameterReader.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include "toml.hpp"

#include <QDebug>
#include <QJsonValue>

#include <iostream>
#include <sstream>



namespace PythiaParameterReaderConstants {
    const QString kMzMinMS2  = QStringLiteral("mzMinMS2");
    const QString kMzMaxMS2  = QStringLiteral("mzMaxMS2");

    const QString kChargeStateMin = QStringLiteral("chargeStateMin");
    const QString kChargeStateMax = QStringLiteral("chargeStateMax");
    const QString kMS2ExtractionWidthPPM = QStringLiteral("ms2ExtractionWidthPPM");
    const QString kPrecursorExtractionWindowThomsons = QStringLiteral("precursorExtractionWindowThomsons");
    const QString kPercentFDR = QStringLiteral("percentFDR");

    const QString kSkipScanCount = QStringLiteral("skipScanCount");
    const QString kMinScanCount = QStringLiteral("minScanCount");
    const QString kUseMeanMz = QStringLiteral("useMeanMz");
    const QString kFilterLength = QStringLiteral("filterLength");
    const QString kSmoothCount = QStringLiteral("smoothCount");
    const QString kSigma = QStringLiteral("sigma");
    const QString kSignalToNoiseRatio = QStringLiteral("signalToNoiseRatio");
    const QString kStopThresholdFraction = QStringLiteral("stopThresholdFraction");
    const QString kScanTimeWindowStDevs = QStringLiteral("scanTimeWindowStDevs");
    const QString kReportDecoys = QStringLiteral("reportDecoys");
    const QString kSubtractShadows = QStringLiteral("subtractShadows");
    const QString kThreadCount = QStringLiteral("threadCount");

    const QString kDigestParams = QStringLiteral("DigestParams");
    const QString kMS2Params = QStringLiteral("MS2Params");
    const QString kPrecursorParams = QStringLiteral("PrecursorParams");
    const QString kPeakIntegrationParams = QStringLiteral("PeakIntegrationParams");
    const QString kFdrParams = QStringLiteral("FdrParams");
    const QString kGeneral = QStringLiteral("General");

}

PythiaParameters PythiaParameterReader::genericPythiaParametersForTests() {

    PythiaParameters pythiaParameters;

    pythiaParameters.ms2ExtractionWidthPPM = 20.0;
    pythiaParameters.precursorExtractionWindowThomsons = 0.0;
    pythiaParameters.chargeStateMin = 2;
    pythiaParameters.chargeStateMax = 3;
    pythiaParameters.minScanCount = 2;
    pythiaParameters.skipScanCount = 0;
    pythiaParameters.filterLength = 5;
    pythiaParameters.smoothCount = 2;
    pythiaParameters.sigma = 1;
    pythiaParameters.signalToNoiseRatio = 2;
    pythiaParameters.mzMinMS2 = 200.0;
    pythiaParameters.mzMaxMS2 = 1500.0;

    return pythiaParameters;
}

Err PythiaParameterReader::buildPythiaParameters(
        const QString &pythiaParametersFilePath,
        PythiaParameters *pythiaParameters
        ) {

    ERR_INIT

    //TODO validate these better.

    using namespace PythiaParameterReaderConstants;

    e = ErrorUtils::fileExists(pythiaParametersFilePath); ree;
    auto parser = toml::parse_file(pythiaParametersFilePath.toStdString());

    const auto ms2ParamsNode =  parser[kMS2Params.toStdString()];
    pythiaParameters->ms2ExtractionWidthPPM = ms2ParamsNode[kMS2ExtractionWidthPPM.toStdString()].value_or(0.0);
    pythiaParameters->mzMinMS2 = ms2ParamsNode[kMzMinMS2.toStdString()].value_or(0.0);
    pythiaParameters->mzMaxMS2 = ms2ParamsNode[kMzMaxMS2.toStdString()].value_or(0.0);

    const auto precursorParamsNode =  parser[kPrecursorParams.toStdString()];
    pythiaParameters->chargeStateMin = precursorParamsNode[kChargeStateMin.toStdString()].value_or(0);
    pythiaParameters->chargeStateMax = precursorParamsNode[kChargeStateMax.toStdString()].value_or(0);
    pythiaParameters->precursorExtractionWindowThomsons = precursorParamsNode[kPrecursorExtractionWindowThomsons.toStdString()].value_or(0.0);

    const auto peakIntegrationParamsNode =  parser[kPeakIntegrationParams.toStdString()];
    pythiaParameters->filterLength = peakIntegrationParamsNode[kFilterLength.toStdString()].value_or(0);
    pythiaParameters->sigma = peakIntegrationParamsNode[kSigma.toStdString()].value_or(0.0);
    pythiaParameters->signalToNoiseRatio = peakIntegrationParamsNode[kSignalToNoiseRatio.toStdString()].value_or(0.0);
    pythiaParameters->smoothCount = peakIntegrationParamsNode[kSmoothCount.toStdString()].value_or(0);

    const auto fdrParamsNode = parser[kFdrParams.toStdString()];
    pythiaParameters->percentFDR = fdrParamsNode[kPercentFDR.toStdString()].value_or(1.0);
    pythiaParameters->reportDecoys = fdrParamsNode[kReportDecoys.toStdString()].value_or(false);

    constexpr int defaultThreadCount = 8;
    const auto generalNode = parser[kGeneral.toStdString()];
    pythiaParameters->threadCount = pythiaParameters->threadCount == 0
                                  ? ParallelUtils::numberOfAvailableSystemProcessors()
                                  : generalNode[kThreadCount.toStdString()].value_or(defaultThreadCount);



    ERR_RETURN
}
