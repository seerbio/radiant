//
// Created by Drucifer on 12/4/2021.
//

#include "PythiaParameterReader.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include "toml.hpp"

#include <QDebug>

#include <iostream>


namespace PythiaParameterReaderConstants {

    const QString kGeneral = QStringLiteral("General");
    const QString kThreadCount = QStringLiteral("threadCount");
    const QString  kVerbosity = QStringLiteral("verbosity");

    const QString kLibraryParams = QStringLiteral("LibraryParams");
    const QString kChargeStateMin = QStringLiteral("chargeStateMin");
    const QString kChargeStateMax = QStringLiteral("chargeStateMax");
    const QString kMzMinMS2  = QStringLiteral("mzMinMS2");
    const QString kMzMaxMS2  = QStringLiteral("mzMaxMS2");
    const QString kPeptideLengthMin = QStringLiteral("peptideLengthMin");
    const QString kPeptideLengthMax  = QStringLiteral("peptideLengthMax");
    const QString kTrancheSizeMax  = QStringLiteral("trancheSizeMax");

    const QString kMS1Params = QStringLiteral("MS1Params");
    const QString kPrecursorExtractionWindowThomsons = QStringLiteral("precursorExtractionWindowThomsons");
    const QString kMS1ExtractionWidthPPM = QStringLiteral("ms1ExtractionWidthPPM");

    const QString kMS2Params = QStringLiteral("MS2Params");
    const QString kFilterLengthIntegration = QStringLiteral("filterLengthIntegration");
    const QString kFilterLengthMS2 = QStringLiteral("filterLengthMS2");
    const QString kIonsSharedToReject = QStringLiteral("ionsSharedToReject");
    const QString kMS2ExtractionWidthPPM = QStringLiteral("ms2ExtractionWidthPPM");
    const QString kMinMs2FragCount = QStringLiteral("minMs2FragCount");
    const QString kScanTimeWindowStDevs = QStringLiteral("scanTimeWindowStDevs");
    const QString kSubtractShadows = QStringLiteral("subtractShadows");
    const QString kSmoothCountMS2 = QStringLiteral("smoothCountMS2");
    const QString kStopThresholdFractionMS2 = QStringLiteral("stopThresholdFractionMS2");

    const QString kFdrParams = QStringLiteral("FdrParams");
    const QString kPercentFDR = QStringLiteral("percentFDR");
    const QString kReportDecoys = QStringLiteral("reportDecoys");

    const QString kPeakIntegrationParams = QStringLiteral("PeakIntegrationParams");
    const QString kFilterLength = QStringLiteral("filterLength");
    const QString kSmoothCount = QStringLiteral("smoothCount");
    const QString kSigma = QStringLiteral("sigma");
    const QString kSignalToNoiseRatio = QStringLiteral("signalToNoiseRatio");
    const QString kStopThresholdFraction = QStringLiteral("stopThresholdFraction");

    const QString kFeatureFinderParams = QStringLiteral("FeatureFinderParams");
    const QString kMinScanCount = QStringLiteral("minScanCount");
    const QString kSkipScanCount = QStringLiteral("skipScanCount");
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

    constexpr int defaultThreadCount = 8;
    const auto generalNode = parser[kGeneral.toStdString()];
    pythiaParameters->threadCount = pythiaParameters->threadCount == 0
                                  ? ParallelUtils::numberOfAvailableSystemProcessors()
                                  : generalNode[kThreadCount.toStdString()].value_or(defaultThreadCount);
    pythiaParameters->verbosity = parser[kGeneral.toStdString()][kVerbosity.toStdString()].value_or(0);

    const auto libraryNode =  parser[kLibraryParams.toStdString()];
    pythiaParameters->chargeStateMin = libraryNode[kChargeStateMin.toStdString()].value_or(0);
    pythiaParameters->chargeStateMax = libraryNode[kChargeStateMax.toStdString()].value_or(0);
    pythiaParameters->mzMinMS2 = libraryNode[kMzMinMS2.toStdString()].value_or(0.0);
    pythiaParameters->mzMaxMS2 = libraryNode[kMzMaxMS2.toStdString()].value_or(0.0);
    pythiaParameters->peptideLengthMin = libraryNode[kPeptideLengthMin.toStdString()].value_or(0);
    pythiaParameters->peptideLengthMax = libraryNode[kPeptideLengthMax.toStdString()].value_or(0);
    pythiaParameters->trancheSizeMax = libraryNode[kTrancheSizeMax.toStdString()].value_or(0);

    const auto ms1ParamsNode =  parser[kMS1Params.toStdString()];
    pythiaParameters->ms1ExtractionWidthPPM = ms1ParamsNode[kMS1ExtractionWidthPPM.toStdString()].value_or(0.0);
    pythiaParameters->precursorExtractionWindowThomsons = ms1ParamsNode[kPrecursorExtractionWindowThomsons.toStdString()].value_or(0.0);

    const auto ms2ParamsNode =  parser[kMS2Params.toStdString()];
    pythiaParameters->filterLengthIntegration = ms2ParamsNode[kFilterLengthIntegration.toStdString()].value_or(0);
    pythiaParameters->filterLengthMS2 = ms2ParamsNode[kFilterLengthMS2.toStdString()].value_or(0);
    pythiaParameters->ionsSharedToReject = ms2ParamsNode[kIonsSharedToReject.toStdString()].value_or(0);
    pythiaParameters->ms2ExtractionWidthPPM = ms2ParamsNode[kMS2ExtractionWidthPPM.toStdString()].value_or(0.0);
    pythiaParameters->minMs2FragCount = ms2ParamsNode[kMinMs2FragCount.toStdString()].value_or(0);
    pythiaParameters->scanTimeWindowStDevs = ms2ParamsNode[kScanTimeWindowStDevs.toStdString()].value_or(0);
    pythiaParameters->subtractShadows = ms2ParamsNode[kSubtractShadows.toStdString()].value_or(true);
    pythiaParameters->smoothCountMS2 = ms2ParamsNode[kSmoothCountMS2.toStdString()].value_or(0);
    pythiaParameters->stopThresholdFractionMS2 = ms2ParamsNode[kStopThresholdFractionMS2.toStdString()].value_or(0.0f);

    const auto fdrParamsNode = parser[kFdrParams.toStdString()];
    pythiaParameters->percentFDR = fdrParamsNode[kPercentFDR.toStdString()].value_or(1.0);
    pythiaParameters->reportDecoys = fdrParamsNode[kReportDecoys.toStdString()].value_or(false);

    const auto peakIntegrationParamsNode =  parser[kPeakIntegrationParams.toStdString()];
    pythiaParameters->filterLength = peakIntegrationParamsNode[kFilterLength.toStdString()].value_or(0);
    pythiaParameters->sigma = peakIntegrationParamsNode[kSigma.toStdString()].value_or(0.0);
    pythiaParameters->signalToNoiseRatio = peakIntegrationParamsNode[kSignalToNoiseRatio.toStdString()].value_or(0.0);
    pythiaParameters->smoothCount = peakIntegrationParamsNode[kSmoothCount.toStdString()].value_or(0);
    pythiaParameters->stopThresholdFraction = ms2ParamsNode[kStopThresholdFraction.toStdString()].value_or(0.0f);

    const auto featureFinderParamsNode = parser[kFeatureFinderParams.toStdString()];
    pythiaParameters->minScanCount = featureFinderParamsNode[kMinScanCount.toStdString()].value_or(0);
    pythiaParameters->skipScanCount = featureFinderParamsNode[kSkipScanCount.toStdString()].value_or(0);

    ERR_RETURN
}
