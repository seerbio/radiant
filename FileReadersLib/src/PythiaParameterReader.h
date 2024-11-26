//
// Created by Drucifer on 12/4/2021.
//

#ifndef PYTHIAPARAMETERREADER_H
#define PYTHIAPARAMETERREADER_H

#include "AminoAcids.h"
#include "Error.h"
#include "FileReadersLib_Exports.h"


using namespace Error;

namespace PythiaParameterReaderConstants {

    extern const QString FILEREADERSLIB_EXPORTS kGeneral;
    extern const QString FILEREADERSLIB_EXPORTS kThreadCount;
    extern const QString FILEREADERSLIB_EXPORTS kVerbosity;
    extern const QString FILEREADERSLIB_EXPORTS kWritePythiaDIA;
    extern const QString FILEREADERSLIB_EXPORTS kUseLazyLoading;

    extern const QString FILEREADERSLIB_EXPORTS kLibraryParams;
    extern const QString FILEREADERSLIB_EXPORTS kChargeStateMin;
    extern const QString FILEREADERSLIB_EXPORTS kChargeStateMax;
    extern const QString FILEREADERSLIB_EXPORTS kMzMinMS2;
    extern const QString FILEREADERSLIB_EXPORTS kMzMaxMS2;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMin;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMax;
    extern const QString FILEREADERSLIB_EXPORTS kTrancheSizeMax;

    extern const QString FILEREADERSLIB_EXPORTS kMS1Params;
    extern const QString FILEREADERSLIB_EXPORTS kPrecursorExtractionWindowThomsons;
    extern const QString FILEREADERSLIB_EXPORTS kMS1ExtractionWidthPPM;
    extern const QString FILEREADERSLIB_EXPORTS kMS1ExtractionWidthPPMOverride;

    extern const QString FILEREADERSLIB_EXPORTS kMS2Params;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLengthIntegration;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLengthMS2;
    extern const QString FILEREADERSLIB_EXPORTS kIonsSharedToReject;
    extern const QString FILEREADERSLIB_EXPORTS kMS2ExtractionWidthPPM;
    extern const QString FILEREADERSLIB_EXPORTS kMS2ExtractionWidthPPMOverride;
    extern const QString FILEREADERSLIB_EXPORTS kMinMs2FragCount;
    extern const QString FILEREADERSLIB_EXPORTS kScanTimeWindowStDevs;
    extern const QString FILEREADERSLIB_EXPORTS kSubtractShadows;
    extern const QString FILEREADERSLIB_EXPORTS kSmoothCountMS2;
    extern const QString FILEREADERSLIB_EXPORTS kStopThresholdFractionMS2;
    extern const QString FILEREADERSLIB_EXPORTS kCalibrationTrainingVolume;
    extern const QString FILEREADERSLIB_EXPORTS kPeakCenter;

    extern const QString FILEREADERSLIB_EXPORTS kTopNIntegrations;
    extern const QString FILEREADERSLIB_EXPORTS kMaxAnchorColumnIndex;

    extern const QString FILEREADERSLIB_EXPORTS kFdrParams;
    extern const QString FILEREADERSLIB_EXPORTS kPercentFDR;
    extern const QString FILEREADERSLIB_EXPORTS kReportDecoys;

    extern const QString FILEREADERSLIB_EXPORTS kPeakIntegrationParams;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLength;
    extern const QString FILEREADERSLIB_EXPORTS kSmoothCount;
    extern const QString FILEREADERSLIB_EXPORTS kSigma;
    extern const QString FILEREADERSLIB_EXPORTS kSignalToNoiseRatio;
    extern const QString FILEREADERSLIB_EXPORTS kStopThresholdFraction;

    extern const QString FILEREADERSLIB_EXPORTS kFeatureFinderParams;
    extern const QString FILEREADERSLIB_EXPORTS kMinScanCount;
    extern const QString FILEREADERSLIB_EXPORTS kSkipScanCount;

    extern const QString FILEREADERSLIB_EXPORTS kRtBinning;


}

struct PythiaParameters{

    //[General]
    int threadCount = 8;
    int verbosity = 0;
    bool writePythiaDIA = true;
    bool useLazyLoading = false;

    //[LibraryParams]
    int chargeStateMin = -1;
    int chargeStateMax = -1;
    double mzMinMS2 = 300.0;
    double mzMaxMS2 = 1999.0;
    int peptideLengthMin = 7;
    int peptideLengthMax = 30;
    int trancheSizeMax = 2e4;

    //[MS1Params]
    double ms1ExtractionWidthPPM = 20.0;
    double ms1ExtractionWidthPPMOverride = -1.0;
    double precursorExtractionWindowThomsons = 0.0;

    //[MS2Params]
    int filterLengthIntegration = 5;
    int filterLengthMS2 = 3;
    int ionsSharedToReject = 4;
    double ms2ExtractionWidthPPM = 20.0;
    double ms2ExtractionWidthPPMOverride = -1.0;
    int minMs2FragCount = 2;
    int rtBinning = 20;
    int scanTimeWindowStDevs = 3;
    bool subtractShadows = true;
    int smoothCountMS2 = 1;
    float stopThresholdFractionMS2 = 0.65;

    int topNIntegrations = 10;
    int maxAnchorColumnIndex = 6;
    int calibrationTrainingVolume = 1000;
    int peakCenter = -1;

    //[FdrParams]
    double percentFDR = 1.0;
    bool reportDecoys = false;

    //[PeakIntegrationParams]
    int filterLength = -1;
    double sigma = -1.0;
    double signalToNoiseRatio = 2.0;
    int smoothCount = -1;
    float stopThresholdFraction = 0.65;

    //[FeatureFinder]
    int minScanCount = 3;
    int skipScanCount = 2;

    [[nodiscard]] bool isValid() const {

        if (chargeStateMin > chargeStateMax) {
            print();
            qDebug()
            << chargeStateMin
            << chargeStateMax
            << PythiaParameterReaderConstants::kChargeStateMin
            << "greater than"
            << PythiaParameterReaderConstants::kChargeStateMax;

            return false;
        }
        if (chargeStateMin < 1) {
            print();
            qDebug() << chargeStateMin << PythiaParameterReaderConstants::kChargeStateMin;
            return false;
        }
        if (chargeStateMax < 1) {
            print();
            qDebug() << chargeStateMax << PythiaParameterReaderConstants::kChargeStateMax;
            return false;
        }
        if (ms2ExtractionWidthPPM < 0) {
            print();
            qDebug() << ms2ExtractionWidthPPM << PythiaParameterReaderConstants::kMS2ExtractionWidthPPM;
            return false;
        }

        if (topNIntegrations <= 0) {
            print();
            qDebug() << topNIntegrations << PythiaParameterReaderConstants::kTopNIntegrations;
            return false;
        }

        if (maxAnchorColumnIndex <= 0) {
            print();
            qDebug() << maxAnchorColumnIndex << PythiaParameterReaderConstants::kMaxAnchorColumnIndex;
            return false;
        }

        return true;
    }

    void print() const {
        qDebug() << QStringLiteral("** Pythia Parameters **************************");

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kGeneral << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kThreadCount) << threadCount;
        qDebug() << qPrintable(PythiaParameterReaderConstants:: kVerbosity) << verbosity;
        qDebug() << qPrintable(PythiaParameterReaderConstants:: kWritePythiaDIA) << writePythiaDIA;
        qDebug() << qPrintable(PythiaParameterReaderConstants:: kUseLazyLoading) << useLazyLoading;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kLibraryParams << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kChargeStateMin) << chargeStateMin;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kChargeStateMax) << chargeStateMax;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMzMinMS2) << mzMinMS2;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMzMaxMS2) << mzMaxMS2;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kPeptideLengthMin) << peptideLengthMin;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kPeptideLengthMax) << peptideLengthMax;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kTrancheSizeMax) << trancheSizeMax;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kMS1Params << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kPrecursorExtractionWindowThomsons) << precursorExtractionWindowThomsons;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMS1ExtractionWidthPPM) << ms1ExtractionWidthPPM;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMS1ExtractionWidthPPMOverride) << ms1ExtractionWidthPPMOverride;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kMS2Params << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kCalibrationTrainingVolume) << calibrationTrainingVolume;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kFilterLengthIntegration) << filterLengthIntegration;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kFilterLengthMS2) << filterLengthMS2;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kIonsSharedToReject) << ionsSharedToReject;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMS2ExtractionWidthPPM) << ms2ExtractionWidthPPM;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMS2ExtractionWidthPPMOverride) << ms2ExtractionWidthPPMOverride;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMinMs2FragCount) << minMs2FragCount;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMaxAnchorColumnIndex) << maxAnchorColumnIndex;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kRtBinning) << rtBinning;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kScanTimeWindowStDevs) << scanTimeWindowStDevs;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSmoothCountMS2) << smoothCountMS2;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kStopThresholdFractionMS2) << stopThresholdFractionMS2;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSubtractShadows) << subtractShadows;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kTopNIntegrations) << topNIntegrations;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kPeakCenter) << peakCenter;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kFdrParams << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kPercentFDR) << percentFDR;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kReportDecoys) << reportDecoys;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kPeakIntegrationParams << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kFilterLength) << filterLength;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSmoothCount) << smoothCount;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSigma) << sigma;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSignalToNoiseRatio) << signalToNoiseRatio;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kStopThresholdFraction) << stopThresholdFraction;

        qDebug() << qPrintable("***") << PythiaParameterReaderConstants::kFeatureFinderParams << qPrintable("***");
        qDebug() << qPrintable(PythiaParameterReaderConstants::kMinScanCount) << minScanCount;
        qDebug() << qPrintable(PythiaParameterReaderConstants::kSkipScanCount) << skipScanCount;

        qDebug() << QStringLiteral("**********************************************");
    }
};


class FILEREADERSLIB_EXPORTS PythiaParameterReader {

public:

    PythiaParameterReader() = default;
    ~PythiaParameterReader() = default;

    static PythiaParameters genericPythiaParametersForTests();

    static Err buildPythiaParameters(
            const QString &pythiaParametersFilePath,
            PythiaParameters *pythiaParameters
    );

};


#endif //PYTHIAPARAMETERREADER_H
