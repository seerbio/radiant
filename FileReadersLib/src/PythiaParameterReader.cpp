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
    const QString kVerbosity = QStringLiteral("verbosity");
    const QString kWritePythiaDIA = QStringLiteral("writePythiaDIA");
    const QString kUseLazyLoading = QStringLiteral("useLazyLoading");
    const QString kReannotate = QStringLiteral("reannotate");
    const QString kShortReport = QStringLiteral("shortReport");

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
    const QString kMS1ExtractionWidthPPMOverride = QStringLiteral("ms1ExtractionWidthPPMOverride");

    const QString kMS2Params = QStringLiteral("MS2Params");
    const QString kFilterLengthIntegration = QStringLiteral("filterLengthIntegration");
    const QString kFilterLengthMS2 = QStringLiteral("filterLengthMS2");
    const QString kIonsSharedToReject = QStringLiteral("ionsSharedToReject");
    const QString kMS2ExtractionWidthPPM = QStringLiteral("ms2ExtractionWidthPPM");
    const QString kMS2ExtractionWidthPPMOverride = QStringLiteral("ms2ExtractionWidthPPMOverride");
    const QString kMinMs2FragCount = QStringLiteral("minMs2FragCount");
    const QString kScanTimeWindowStDevs = QStringLiteral("scanTimeWindowStDevs");
    const QString kSubtractShadows = QStringLiteral("subtractShadows");
    const QString kSmoothCountMS2 = QStringLiteral("smoothCountMS2");
    const QString kStopThresholdFractionMS2 = QStringLiteral("stopThresholdFractionMS2");
    const QString kPeakDifferenceThresholdFraction = QStringLiteral("peakDifferenceThresholdFraction");
    const QString kCalibrationTrainingVolume = QStringLiteral("calibrationTrainingVolume");
    const QString kPeakCenter = QStringLiteral("peakCenter");

    const QString kTopNIntegrations = QStringLiteral("topNIntegrations");
    const QString kMaxAnchorColumnIndex = QStringLiteral("maxAnchorColumnIndex");

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

    const QString kRtBinning = QStringLiteral("rtBinning");

    const QString kNeuralNetParams = QStringLiteral("NeuralNetParams");
    const QString kEpochs = QStringLiteral("epochs");
    const QString kBaggingSize = QStringLiteral("baggingSize");
    const QString kLearningRate = QStringLiteral("learningRate");
    const QString kNodesFraction = QStringLiteral("nodesFraction");
    const QString kFocalLossGamma = QStringLiteral("focalLossGamma");
    const QString kParallelNeuralNets = QStringLiteral("parallelNeuralNets");

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

    const auto generalNode = parser[kGeneral.toStdString()];
    pythiaParameters->threadCount = pythiaParameters->threadCount == 0
                                  ? ParallelUtils::numberOfAvailableSystemProcessors()
                                  : generalNode[kThreadCount.toStdString()].value_or(ParallelUtils::numberOfAvailableSystemProcessors());
    pythiaParameters->verbosity = parser[kGeneral.toStdString()][kVerbosity.toStdString()].value_or(0);
    pythiaParameters->useLazyLoading = parser[kGeneral.toStdString()][kUseLazyLoading.toStdString()].value_or(false);
    pythiaParameters->writePythiaDIA = parser[kGeneral.toStdString()][kWritePythiaDIA.toStdString()].value_or(true);
    pythiaParameters->reannotate = parser[kGeneral.toStdString()][kReannotate.toStdString()].value_or(false);
    pythiaParameters->shortReport = parser[kGeneral.toStdString()][kShortReport.toStdString()].value_or(false);

    const auto libraryNode =  parser[kLibraryParams.toStdString()];
    pythiaParameters->chargeStateMin = libraryNode[kChargeStateMin.toStdString()].value_or(0);
    pythiaParameters->chargeStateMax = libraryNode[kChargeStateMax.toStdString()].value_or(0);
    pythiaParameters->mzMinMS2 = libraryNode[kMzMinMS2.toStdString()].value_or(0.0);
    pythiaParameters->mzMaxMS2 = libraryNode[kMzMaxMS2.toStdString()].value_or(0.0);
    pythiaParameters->peptideLengthMin = libraryNode[kPeptideLengthMin.toStdString()].value_or(0);
    pythiaParameters->peptideLengthMax = libraryNode[kPeptideLengthMax.toStdString()].value_or(0);
    pythiaParameters->trancheSizeMax = libraryNode[kTrancheSizeMax.toStdString()].value_or(0);

    const auto ms1ParamsNode =  parser[kMS1Params.toStdString()];
    pythiaParameters->ms1ExtractionWidthPPM = ms1ParamsNode[kMS1ExtractionWidthPPM.toStdString()].value_or(pythiaParameters->ms1ExtractionWidthPPM);
    pythiaParameters->ms1ExtractionWidthPPMOverride = ms1ParamsNode[kMS1ExtractionWidthPPMOverride.toStdString()].value_or(-1.0);
    pythiaParameters->precursorExtractionWindowThomsons = ms1ParamsNode[kPrecursorExtractionWindowThomsons.toStdString()].value_or(pythiaParameters->precursorExtractionWindowThomsons);

    const auto ms2ParamsNode =  parser[kMS2Params.toStdString()];
    pythiaParameters->filterLengthIntegration = ms2ParamsNode[kFilterLengthIntegration.toStdString()].value_or(pythiaParameters->filterLengthIntegration);
    pythiaParameters->filterLengthMS2 = ms2ParamsNode[kFilterLengthMS2.toStdString()].value_or(pythiaParameters->filterLengthMS2);
    pythiaParameters->ionsSharedToReject = ms2ParamsNode[kIonsSharedToReject.toStdString()].value_or(pythiaParameters->ionsSharedToReject);
    pythiaParameters->ms2ExtractionWidthPPM = ms2ParamsNode[kMS2ExtractionWidthPPM.toStdString()].value_or(pythiaParameters->ms2ExtractionWidthPPM);
    pythiaParameters->ms2ExtractionWidthPPMOverride = ms2ParamsNode[kMS2ExtractionWidthPPMOverride.toStdString()].value_or(-1.0);
    pythiaParameters->minMs2FragCount = ms2ParamsNode[kMinMs2FragCount.toStdString()].value_or(pythiaParameters->minMs2FragCount);
    pythiaParameters->scanTimeWindowStDevs = ms2ParamsNode[kScanTimeWindowStDevs.toStdString()].value_or(pythiaParameters->scanTimeWindowStDevs);
    pythiaParameters->subtractShadows = ms2ParamsNode[kSubtractShadows.toStdString()].value_or(true);
    pythiaParameters->smoothCountMS2 = ms2ParamsNode[kSmoothCountMS2.toStdString()].value_or(pythiaParameters->smoothCountMS2);
    pythiaParameters->stopThresholdFractionMS2 = ms2ParamsNode[kStopThresholdFractionMS2.toStdString()].value_or(pythiaParameters->stopThresholdFractionMS2);
    pythiaParameters->peakDifferenceFractionThreshold = ms2ParamsNode[kPeakDifferenceThresholdFraction.toStdString()].value_or(pythiaParameters->peakDifferenceFractionThreshold);
    pythiaParameters->calibrationTrainingVolume = ms2ParamsNode[kCalibrationTrainingVolume.toStdString()].value_or(pythiaParameters->calibrationTrainingVolume);
    pythiaParameters->rtBinning = ms2ParamsNode[kRtBinning.toStdString()].value_or(pythiaParameters->rtBinning);
    pythiaParameters->peakCenter = ms2ParamsNode[kPeakCenter.toStdString()].value_or(pythiaParameters->peakCenter);
    pythiaParameters->topNIntegrations = ms2ParamsNode[kTopNIntegrations.toStdString()].value_or(pythiaParameters->topNIntegrations);
    pythiaParameters->maxAnchorColumnIndex = ms2ParamsNode[kMaxAnchorColumnIndex.toStdString()].value_or(12);

    const auto fdrParamsNode = parser[kFdrParams.toStdString()];
    pythiaParameters->percentFDR = fdrParamsNode[kPercentFDR.toStdString()].value_or(1.0);
    pythiaParameters->reportDecoys = fdrParamsNode[kReportDecoys.toStdString()].value_or(false);

    const auto peakIntegrationParamsNode =  parser[kPeakIntegrationParams.toStdString()];
    pythiaParameters->filterLength = peakIntegrationParamsNode[kFilterLength.toStdString()].value_or(0);
    pythiaParameters->sigma = peakIntegrationParamsNode[kSigma.toStdString()].value_or(0.0);
    pythiaParameters->signalToNoiseRatio = peakIntegrationParamsNode[kSignalToNoiseRatio.toStdString()].value_or(0.0);
    pythiaParameters->smoothCount = peakIntegrationParamsNode[kSmoothCount.toStdString()].value_or(0);
    pythiaParameters->stopThresholdFraction = ms2ParamsNode[kStopThresholdFraction.toStdString()].value_or(pythiaParameters->stopThresholdFraction);

    const auto featureFinderParamsNode = parser[kFeatureFinderParams.toStdString()];
    pythiaParameters->minScanCount = featureFinderParamsNode[kMinScanCount.toStdString()].value_or(0);
    pythiaParameters->skipScanCount = featureFinderParamsNode[kSkipScanCount.toStdString()].value_or(0);

    const auto neuralNetParamsNode = parser[kNeuralNetParams.toStdString()];
    pythiaParameters->epochs = neuralNetParamsNode[kEpochs.toStdString()].value_or(pythiaParameters->epochs);
    pythiaParameters->baggingSize = neuralNetParamsNode[kBaggingSize.toStdString()].value_or(pythiaParameters->baggingSize);
    pythiaParameters->learningRate = neuralNetParamsNode[kLearningRate.toStdString()].value_or(pythiaParameters->learningRate);
    pythiaParameters->nodesFraction = neuralNetParamsNode[kNodesFraction.toStdString()].value_or(pythiaParameters->nodesFraction);
    pythiaParameters->focalLossGamma = neuralNetParamsNode[kFocalLossGamma.toStdString()].value_or(pythiaParameters->focalLossGamma);
    pythiaParameters->parallelNeuralNets = neuralNetParamsNode[kParallelNeuralNets.toStdString()].value_or(pythiaParameters->parallelNeuralNets);

	if (pythiaParameters->baggingSize < 2) {
		qDebug()
		<< qPrintable(S_GLOBAL_TIMER.elapsed())
		<< "Bagging size must be at least 2.  Bagging size has been adjusted from"
		<< pythiaParameters->baggingSize
		<< "to 2.";
		pythiaParameters->baggingSize = 2;
	}

    ERR_RETURN
}
