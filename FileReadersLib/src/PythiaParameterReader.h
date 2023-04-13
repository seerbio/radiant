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
    extern const QString FILEREADERSLIB_EXPORTS kAddDecoys;
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

    int positionalLocationIndexes(const QString &positionalLocation,
                                  const QString &peptideSequence);

    static QStringList modKeys();

    void print() const;
};


struct PythiaParameters{
    //TODO make an is valid method to check and see if all values are initiated.
    int trancheSize = 16;

    QStringList nTermCleavePoints;
    QStringList cTermCleavePoints;

    Raggedness raggedness = Raggedness::NoRagged;

    int allowedMissedCleavages = 0;
    int peptideLengthMin = 7;
    int peptideLengthMax = 50;
    int maxModificationsPeptide = 1;

    int chargeStateMin = -1;
    int chargeStateMax = -1;
    int maxTandemPointCount = -1;
    int returnPSMTopN = -1;

    double ms2ExtractionWidthPPM = -1.0;
    double precursorExtractionWindowThomsons = -1.0;
    double percentFDR = 1.0;
    double mzMinDataStructure = 300.0;
    double mzMaxDataStructure = 1999.0;

    bool replaceLeucinesWithX = true;
    bool addDecoys = true; //TODO change this to a number for rounds of decoys.

    //TODO hook these up
    double featureFinderTolerancePPM = 12;
    int skipScanCount = 2;
    int minScanCount = 3;
    bool useMeanMz = true;

    //TODO hook these up
    bool denoise = true;
    bool deisotope = true;
    bool smooth = false;

    //TODO hook these up
    int topNMs2Ions = 12;
    int minFoundMzPeaks = 5;

    [[nodiscard]] bool isValid() const {

        if (chargeStateMin > chargeStateMax) {
            return false;
        }
        if (chargeStateMin < 1) {
            return false;
        }
        if (chargeStateMax < 1) {
            return false;
        }
        if (maxTandemPointCount < 0) {
            return false;
        }
        if (returnPSMTopN < 0) {
            return false;
        }
        if (ms2ExtractionWidthPPM < 0) {
            return false;
        }
        if (precursorExtractionWindowThomsons < 0) {
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
        qDebug() << PythiaParameterReaderConstants::kAddDecoys << addDecoys;
        qDebug() << "topNMs2Ions" << topNMs2Ions;
        qDebug() << featureFinderTolerancePPM << skipScanCount << minScanCount << useMeanMz; //TODO make this proper like the rest

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

    static Err applyFixedModificationsToAminoAcids(const PythiaParameters &reader,
                                                   AminoAcids *aminoAcids);

private:

    Err validateJsonKeys();

};


#endif //PYTHIAPARAMETERREADER_H
