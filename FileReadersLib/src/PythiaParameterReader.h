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

    extern const QString FILEREADERSLIB_EXPORTS kMS2Params;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLengthIntegration;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLengthMS2;
    extern const QString FILEREADERSLIB_EXPORTS kIonsSharedToReject;
    extern const QString FILEREADERSLIB_EXPORTS kMS2ExtractionWidthPPM;
    extern const QString FILEREADERSLIB_EXPORTS kMinMs2FragCount;
    extern const QString FILEREADERSLIB_EXPORTS kScanTimeWindowStDevs;
    extern const QString FILEREADERSLIB_EXPORTS kSubtractShadows;
    extern const QString FILEREADERSLIB_EXPORTS kSmoothCountMS2;
    extern const QString FILEREADERSLIB_EXPORTS kStopThresholdFractionMS2;

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

}

struct PythiaParameters{

    //[General]
    int threadCount = 8;
    int verbosity = 0;

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
    double precursorExtractionWindowThomsons = 0.0;

    //[MS2Params]
    int filterLengthIntegration = 5;
    int filterLengthMS2 = 3;
    int ionsSharedToReject = 4;
    double ms2ExtractionWidthPPM = 20.0;
    int minMs2FragCount = 2;
    int scanTimeWindowStDevs = 3;
    bool subtractShadows = true;
    int smoothCountMS2 = 1;
    float stopThresholdFractionMS2 = 0.65;

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
            qDebug() << chargeStateMin << chargeStateMax << "Charge state min greater than max";
            return false;
        }
        if (chargeStateMin < 1) {
            print();
            qDebug() << chargeStateMin << "charge state min";
            return false;
        }
        if (chargeStateMax < 1) {
            print();
            qDebug() << chargeStateMax << "charge state max";
            return false;
        }
        if (ms2ExtractionWidthPPM < 0) {
            print();
            qDebug() << ms2ExtractionWidthPPM << "ms2ExtractPPM";
            return false;
        }

        return true;
    }

    void print() const {
        qDebug() << QStringLiteral("** Digest Parameters **************************");

        qDebug() << PythiaParameterReaderConstants::kGeneral;
        qDebug() << PythiaParameterReaderConstants::kThreadCount << threadCount;
        qDebug() << PythiaParameterReaderConstants:: kVerbosity << verbosity;

        qDebug() << PythiaParameterReaderConstants::kLibraryParams;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMin << chargeStateMin;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMax << chargeStateMax;
        qDebug() << PythiaParameterReaderConstants::kMzMinMS2 << mzMinMS2;
        qDebug() << PythiaParameterReaderConstants::kMzMaxMS2 << mzMaxMS2;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMin << peptideLengthMin;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMax << peptideLengthMax;
        qDebug() << PythiaParameterReaderConstants::kTrancheSizeMax << trancheSizeMax;

        qDebug() << PythiaParameterReaderConstants::kMS1Params;
        qDebug() << PythiaParameterReaderConstants::kPrecursorExtractionWindowThomsons << precursorExtractionWindowThomsons;
        qDebug() << PythiaParameterReaderConstants::kMS1ExtractionWidthPPM << ms1ExtractionWidthPPM;

        qDebug() << PythiaParameterReaderConstants::kMS2Params;
        qDebug() << PythiaParameterReaderConstants::kFilterLengthIntegration << filterLengthIntegration;
        qDebug() << PythiaParameterReaderConstants::kFilterLengthMS2 << filterLengthMS2;
        qDebug() << PythiaParameterReaderConstants::kIonsSharedToReject << ionsSharedToReject;
        qDebug() << PythiaParameterReaderConstants::kMS2ExtractionWidthPPM << ms2ExtractionWidthPPM;
        qDebug() << PythiaParameterReaderConstants::kMinMs2FragCount << minMs2FragCount;
        qDebug() << PythiaParameterReaderConstants::kScanTimeWindowStDevs << scanTimeWindowStDevs;
        qDebug() << PythiaParameterReaderConstants::kSubtractShadows << subtractShadows;
        qDebug() << PythiaParameterReaderConstants::kSmoothCountMS2 << smoothCountMS2;
        qDebug() << PythiaParameterReaderConstants::kStopThresholdFractionMS2 << stopThresholdFractionMS2;

        qDebug() << PythiaParameterReaderConstants::kFdrParams;
        qDebug() << PythiaParameterReaderConstants::kPercentFDR << percentFDR;
        qDebug() << PythiaParameterReaderConstants::kReportDecoys << reportDecoys;

        qDebug() << PythiaParameterReaderConstants::kPeakIntegrationParams;
        qDebug() << PythiaParameterReaderConstants::kFilterLength << filterLength;
        qDebug() << PythiaParameterReaderConstants::kSmoothCount << smoothCount;
        qDebug() << PythiaParameterReaderConstants::kSigma << sigma;
        qDebug() << PythiaParameterReaderConstants::kSignalToNoiseRatio << signalToNoiseRatio;
        qDebug() << PythiaParameterReaderConstants::kStopThresholdFraction << stopThresholdFraction;

        qDebug() << PythiaParameterReaderConstants::kFeatureFinderParams;
        qDebug() << PythiaParameterReaderConstants::kMinScanCount << minScanCount;
        qDebug() << PythiaParameterReaderConstants::kSkipScanCount << skipScanCount;

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
