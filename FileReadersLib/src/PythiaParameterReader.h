//
// Created by Drucifer on 12/4/2021.
//

#ifndef PYTHIAPARAMETERREADER_H
#define PYTHIAPARAMETERREADER_H

#include "AminoAcids.h"
#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "JsonParametersReader.h"


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
    extern const QString FILEREADERSLIB_EXPORTS kMzMinDataStructure;
    extern const QString FILEREADERSLIB_EXPORTS kMzMaxDataStructure;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMin;
    extern const QString FILEREADERSLIB_EXPORTS kPeptideLengthMax;
    extern const QString FILEREADERSLIB_EXPORTS kNoRagged;
    extern const QString FILEREADERSLIB_EXPORTS kCTermRagged;
    extern const QString FILEREADERSLIB_EXPORTS kNTermRagged;
    extern const QString FILEREADERSLIB_EXPORTS kBothRagged;

    extern const QString FILEREADERSLIB_EXPORTS kMaxTandemPointCount;
    extern const QString FILEREADERSLIB_EXPORTS kReturnPSMTopN;
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
    extern const QString FILEREADERSLIB_EXPORTS kTopNMs2Ions;
    extern const QString FILEREADERSLIB_EXPORTS kMinFoundMzPeaks;

    extern const QString FILEREADERSLIB_EXPORTS kDeisotopeScans;
    extern const QString FILEREADERSLIB_EXPORTS kTrancheSizeMax;
    extern const QString FILEREADERSLIB_EXPORTS kCosineSimToAnchorThreshold;
    extern const QString FILEREADERSLIB_EXPORTS kScanTimeWindowMinutes;
    extern const QString FILEREADERSLIB_EXPORTS kReportDecoys;

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

    Modification(QChar residue,
                 const QString &name,
                 const ModificationType &type,
                 const QString &formula);

    Modification(const QString &positionalLocation,
                 const QString &name,
                 const ModificationType &type,
                 const QString &formula);

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
    int maxTandemPointCount = -1;
    int returnPSMTopN = -1;

    double precursorExtractionWindowThomsons = 0.0;
    double percentFDR = 1.0;
    double mzMinDataStructure = 300.0;
    double mzMaxDataStructure = 1999.0;

    double fragIntensityThreshold = 0.025;
    double pValThreshold = 0.05;

    int skipScanCount = 2;
    int minScanCount = 3;
    bool useMeanMz = true;
    int filterLength = -1;
    int smoothCount = -1;
    double sigma = -1.0;
    double signalToNoiseRatio = -1.0;
    int topNMs2Ions = -1;
    int minFoundMzPeaks = -1;
    float stopThresholdFraction = 0.2;

    double cosineSimToAnchorThreshold = 0.4;
    bool subtractShadows = true;

    double scanTimeWindowMinutes = -1.0;
    double ms2ExtractionWidthPPM = -1.0;
    int trancheSizeMax = 5e4;

    int verbosity = 1;
    bool reportDecoys = false;

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
        if (maxTandemPointCount < 0) {
            print();
            qDebug() << maxTandemPointCount << "max tandem point count";
            return false;
        }
        if (returnPSMTopN < 0) {
            print();
            qDebug() << returnPSMTopN << "return PSM Top N";
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
        qDebug() << PythiaParameterReaderConstants::kMzMinDataStructure << mzMinDataStructure;
        qDebug() << PythiaParameterReaderConstants::kMzMaxDataStructure << mzMaxDataStructure;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMin << peptideLengthMin;
        qDebug() << PythiaParameterReaderConstants::kPeptideLengthMax << peptideLengthMax;
        qDebug() << PythiaParameterReaderConstants::kMaxTandemPointCount << maxTandemPointCount;
        qDebug() << PythiaParameterReaderConstants::kReturnPSMTopN << returnPSMTopN;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMin << chargeStateMin;
        qDebug() << PythiaParameterReaderConstants::kChargeStateMax << chargeStateMax;
        qDebug() << PythiaParameterReaderConstants::kMS2ExtractionWidthPPM << ms2ExtractionWidthPPM;
        qDebug() << PythiaParameterReaderConstants::kPrecursorExtractionWindowThomsons << precursorExtractionWindowThomsons;
        qDebug() << PythiaParameterReaderConstants::kPercentFDR << percentFDR;
        qDebug() << PythiaParameterReaderConstants::kMaxModificationsPeptide << maxModificationsPeptide;

        qDebug() << PythiaParameterReaderConstants::kSkipScanCount << skipScanCount;
        qDebug() << PythiaParameterReaderConstants::kMinScanCount << minScanCount;
        qDebug() << PythiaParameterReaderConstants::kUseMeanMz << useMeanMz;
        qDebug() << PythiaParameterReaderConstants::kFilterLength << filterLength;
        qDebug() << PythiaParameterReaderConstants::kSmoothCount << smoothCount;
        qDebug() << PythiaParameterReaderConstants::kSigma << sigma;
        qDebug() << PythiaParameterReaderConstants::kSignalToNoiseRatio << signalToNoiseRatio;
        qDebug() << PythiaParameterReaderConstants::kTopNMs2Ions << topNMs2Ions;
        qDebug() << PythiaParameterReaderConstants::kMinFoundMzPeaks << minFoundMzPeaks;
        qDebug() << PythiaParameterReaderConstants::kCosineSimToAnchorThreshold << cosineSimToAnchorThreshold;
        qDebug() << PythiaParameterReaderConstants::kScanTimeWindowMinutes << scanTimeWindowMinutes;
        qDebug() << PythiaParameterReaderConstants::kTrancheSizeMax << trancheSizeMax;
        qDebug() << PythiaParameterReaderConstants::kReportDecoys << reportDecoys;

        qDebug() << PythiaParameterReaderConstants::kModifications;
        for (const Modification &mod : modifications) {
            mod.print();
        }

        qDebug() << QStringLiteral("**********************************************");
    }

};


class FILEREADERSLIB_EXPORTS PythiaParameterReader : public JsonParametersReader {

public:

    PythiaParameterReader() = default;
    ~PythiaParameterReader() = default;

    Err loadPythiaParameters(PythiaParameters *pythiaParameters);

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

    Err validateJsonKeys();

};


#endif //PYTHIAPARAMETERREADER_H
