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
    extern const QString FILEREADERSLIB_EXPORTS kNTermCleavePoints;
    extern const QString FILEREADERSLIB_EXPORTS kCTermCleavePoints;
    extern const QString FILEREADERSLIB_EXPORTS kRaggedness;
    extern const QString FILEREADERSLIB_EXPORTS kAllowedMissedCleavages;
    extern const QString FILEREADERSLIB_EXPORTS kDynamic;
    extern const QString FILEREADERSLIB_EXPORTS kFixed;
    extern const QString FILEREADERSLIB_EXPORTS kMaxModificationsPeptide;
    extern const QString FILEREADERSLIB_EXPORTS kModifications;
    extern const QString FILEREADERSLIB_EXPORTS kMzMinMS2;
    extern const QString FILEREADERSLIB_EXPORTS kMzMaxMS2;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMin;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMax;
    extern const QString FILEREADERSLIB_EXPORTS kNoRagged;
    extern const QString FILEREADERSLIB_EXPORTS kCTermRagged;
    extern const QString FILEREADERSLIB_EXPORTS kNTermRagged;
    extern const QString FILEREADERSLIB_EXPORTS kBothRagged;

    extern const QString FILEREADERSLIB_EXPORTS kChargeStateMin;
    extern const QString FILEREADERSLIB_EXPORTS kChargeStateMax;
    extern const QString FILEREADERSLIB_EXPORTS kMS2ExtractionWidthPPM;
    extern const QString FILEREADERSLIB_EXPORTS kPrecursorExtractionWindowThomsons;
    extern const QString FILEREADERSLIB_EXPORTS kPercentFDR;

    extern const QString FILEREADERSLIB_EXPORTS kSkipScanCount;
    extern const QString FILEREADERSLIB_EXPORTS kMinScanCount;
    extern const QString FILEREADERSLIB_EXPORTS kUseMeanMz;
    extern const QString FILEREADERSLIB_EXPORTS kFilterLength;
    extern const QString FILEREADERSLIB_EXPORTS kSmoothCount;
    extern const QString FILEREADERSLIB_EXPORTS kSigma;
    extern const QString FILEREADERSLIB_EXPORTS kSignalToNoiseRatio;
    extern const QString FILEREADERSLIB_EXPORTS kMinFoundMzPeaks;
    extern const QString FILEREADERSLIB_EXPORTS kStopThresholdFraction;
    extern const QString FILEREADERSLIB_EXPORTS kCosineSimToAnchorThreshold;
    extern const QString FILEREADERSLIB_EXPORTS kScanTimeWindowStDevs;
    extern const QString FILEREADERSLIB_EXPORTS kReportDecoys;

    extern const QString FILEREADERSLIB_EXPORTS kSubtractShadows;

    extern const QString FILEREADERSLIB_EXPORTS kDigestParams;
    extern const QString FILEREADERSLIB_EXPORTS kMS2Params;
    extern const QString FILEREADERSLIB_EXPORTS kPrecursorParams;
    extern const QString FILEREADERSLIB_EXPORTS kPeakIntegrationParams;
    extern const QString FILEREADERSLIB_EXPORTS kFdrParams;

    extern const QString FILEREADERSLIB_EXPORTS kBypassNeuralNet;

}


enum Raggedness {
    NoRagged = 0
    , CTermRagged
    , NTermRagged
    , BothRagged
};


enum class ModificationType {
    DYNAMIC,
    FIXED
};


class FILEREADERSLIB_EXPORTS Modification {

public:

    Modification(
            QChar residue,
            const QString &name,
            const ModificationType &type,
            const QString &formula
            );

    Modification(
            const QString &positionalLocation,
            const QString &name,
            const ModificationType &type,
            const QString &formula
            );

    Modification() = default;
    ~Modification() = default;

    QChar residue;
    QString positionalLocation;
    QString name;
    ModificationType type = ModificationType::DYNAMIC;
    QString formula;

    static QString nTermPeptide();
    static QString nTermProtein();
    static QString cTermPeptide();
    static QString cTermProtein();

    int positionalLocationIndexes(
            const QString &positionalLocation,
            const QString &peptideSequence
            );

    static QStringList modKeys();

    void print() const;
};

struct PythiaParameters{

    QStringList nTermCleavePoints;
    QStringList cTermCleavePoints;

    Raggedness raggedness = Raggedness::NoRagged;

    int allowedMissedCleavages = 0;
    int peptideLengthMin = 7;
    int peptideLengthMax = 35;
    int maxModificationsPeptide = 1;

    int chargeStateMin = -1;
    int chargeStateMax = -1;

    double precursorExtractionWindowThomsons = 0.0;
    double percentFDR = 1.0;
    double mzMinMS2 = 300.0;
    double mzMaxMS2 = 1999.0;

    int skipScanCount = 2;
    int minScanCount = 3;
    double signalToNoiseRatio = 2.0;

    int filterLength = -1;
    int smoothCount = -1;
    double sigma = -1.0;
    float stopThresholdFraction = 0.2;
    int minFoundMzPeaks = -1;

    double cosineSimToAnchorThreshold = 0.4;

    double ms2ExtractionWidthPPM = -1.0;

    int verbosity = 1;
    bool reportDecoys = false;
    bool subtractShadows = true;

    int scanTimeWindowStDevs = 3;
    int trancheSizeMax = 1e4;
    int ionsSharedToReject = 4;
    int topNMs2Ions = 12;

    bool bypassNeuralNet = false;

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
        if (precursorExtractionWindowThomsons < 0) {
            print();
            qDebug() << precursorExtractionWindowThomsons << "precursorExtractWindowThompsons";
            return false;
        }
        if (topNMs2Ions < 6) {
            print();
            qDebug() << topNMs2Ions << "topNMS2Ions";
            return false;
        }
        if (minFoundMzPeaks < 0) {
            print();
            qDebug() << minFoundMzPeaks << "minFoundMzPeaks";
            return false;
        }

        return true;
    }

    QVector<Modification> modifications;

    AminoAcids aminoAcids;

    void print() const {
        qDebug() << QStringLiteral("** Digest Parameters **************************");
        qDebug() << PythiaParameterReaderConstants::kNTermCleavePoints << nTermCleavePoints;
        qDebug() << PythiaParameterReaderConstants::kCTermCleavePoints << cTermCleavePoints;
        qDebug() << PythiaParameterReaderConstants::kRaggedness << static_cast<std::underlying_type<Raggedness>::type>(raggedness);
        qDebug() << PythiaParameterReaderConstants::kAllowedMissedCleavages << allowedMissedCleavages;
        qDebug() << PythiaParameterReaderConstants::kMzMinMS2 << mzMinMS2;
        qDebug() << PythiaParameterReaderConstants::kMzMaxMS2 << mzMaxMS2;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMin << peptideLengthMin;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMax << peptideLengthMax;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMin << chargeStateMin;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMax << chargeStateMax;
        qDebug() << PythiaParameterReaderConstants::kMS2ExtractionWidthPPM << ms2ExtractionWidthPPM;
        qDebug() << PythiaParameterReaderConstants::kPrecursorExtractionWindowThomsons << precursorExtractionWindowThomsons;
        qDebug() << PythiaParameterReaderConstants::kPercentFDR << percentFDR;
        qDebug() << PythiaParameterReaderConstants::kMaxModificationsPeptide << maxModificationsPeptide;

        qDebug() << PythiaParameterReaderConstants::kSkipScanCount << skipScanCount;
        qDebug() << PythiaParameterReaderConstants::kMinScanCount << minScanCount;

        qDebug() << PythiaParameterReaderConstants::kFilterLength << filterLength;
        qDebug() << PythiaParameterReaderConstants::kSmoothCount << smoothCount;
        qDebug() << PythiaParameterReaderConstants::kSigma << sigma;
        qDebug() << PythiaParameterReaderConstants::kSignalToNoiseRatio << signalToNoiseRatio;
        qDebug() << PythiaParameterReaderConstants::kMinFoundMzPeaks << minFoundMzPeaks;
        qDebug() << PythiaParameterReaderConstants::kStopThresholdFraction << stopThresholdFraction;
        qDebug() << PythiaParameterReaderConstants::kCosineSimToAnchorThreshold << cosineSimToAnchorThreshold;
        qDebug() << PythiaParameterReaderConstants::kScanTimeWindowStDevs << scanTimeWindowStDevs;
        qDebug() << PythiaParameterReaderConstants::kReportDecoys << reportDecoys;
        qDebug() << PythiaParameterReaderConstants::kSubtractShadows << subtractShadows;

        qDebug() << PythiaParameterReaderConstants::kBypassNeuralNet << bypassNeuralNet;

        qDebug() << PythiaParameterReaderConstants::kModifications;
        for (const Modification &mod : modifications) {
            mod.print();
        }

        qDebug() << QStringLiteral("**********************************************");
    }

};


class FILEREADERSLIB_EXPORTS PythiaParameterReader {

public:

    PythiaParameterReader() = default;
    ~PythiaParameterReader() = default;

    static Err applyFixedModificationsToAminoAcids(
            const PythiaParameters &reader,
            AminoAcids *aminoAcids
            );

    static PythiaParameters genericPythiaParametersForTests();

    static Err buildPythiaParameters(
            const QString &pythiaParametersFilePath,
            PythiaParameters *pythiaParameters
    );

private:

    Err validateJsonKeys(); //TODO update this for toml

};


#endif //PYTHIAPARAMETERREADER_H
