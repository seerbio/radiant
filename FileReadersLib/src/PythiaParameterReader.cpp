//
// Created by Drucifer on 12/4/2021.
//

#include "PythiaParameterReader.h"
#include "ErrorUtils.h"

#include <QDebug>
#include <QJsonValue>

Modification::Modification(QChar residue,
                           const QString &name,
                           const ModificationType &type,
                           const QString &formula)
: residue(residue)
, name(name)
, type(type)
, formula(formula){}


Modification::Modification(
        const QString &positionalLocation,
        const QString &name,
        const ModificationType &type,
        const QString &formula
        )
: positionalLocation(positionalLocation)
, name(name)
, type(type)
, formula(formula) {}


void Modification::print() const {

    using namespace PythiaParameterReaderConstants;
    const QString modType = type == ModificationType::FIXED ? kFixed : kDynamic;

    const QString locationToPrint = positionalLocation.isEmpty() ? residue : positionalLocation;

    qDebug() << "Residue:" << locationToPrint << "Name:" << name
             << "Mod Type:" << modType << "Formula:" << formula;
}


QStringList Modification::modKeys() {
    return {"residue", "name", "type", "formula"};
}


int Modification::positionalLocationIndexes(
        const QString &positionalLocation,
        const QString &peptideSequence
        ) {


    if (positionalLocation == nTermProtein() || positionalLocation == nTermPeptide()) {
        return 0;
    }

    else if (positionalLocation == cTermProtein() || positionalLocation == cTermPeptide()) {
        return peptideSequence.size() - 1;
    }

    return -1;
}


QString Modification::nTermPeptide() {
    return QStringLiteral("N-term-peptide");
}


QString Modification::nTermProtein() {
    return QStringLiteral("N-term-protein");
}


QString Modification::cTermPeptide() {
    return QStringLiteral("C-term-peptide");
}


QString Modification::cTermProtein() {
    return QStringLiteral("C-term-protein");
}


namespace PythiaParameterReaderConstants {
    const QString kAddDecoys = QStringLiteral("addDecoys");
    const QString kNTermCleavePoints = QStringLiteral("nTermCleavePoints");
    const QString kCTermCleavePoints = QStringLiteral("cTermCleavePoints");
    const QString kRaggedness = QStringLiteral("raggedness");
    const QString kAllowedMissedCleavages = QStringLiteral("allowedMissedCleavages");
    const QString kFixed = QStringLiteral("fixed");
    const QString kDynamic = QStringLiteral("dynamic");
    const QString kMaxModificationsPeptide = QStringLiteral("maxModificationsPeptide");
    const QString kModifications = QStringLiteral("modifications");
    const QString kMzMinDataStructure  = QStringLiteral("mzMinDataStructure");
    const QString kMzMaxDataStructure  = QStringLiteral("mzMaxDataStructure");
    const QString kPeptideLengthMin = QStringLiteral("peptideLengthMin");
    const QString kPeptideLengthMax = QStringLiteral("peptideLengthMax");
    const QString kNoRagged = QStringLiteral("NoRagged");
    const QString kCTermRagged = QStringLiteral("CTermRagged");
    const QString kNTermRagged = QStringLiteral("NTermRagged");
    const QString kBothRagged = QStringLiteral("BothRagged");

    const QString kMaxTandemPointCount  = QStringLiteral("maxTandemPointCount");
    const QString kReturnPSMTopN = QStringLiteral("returnPSMTopN");
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
    const QString kTopNMs2Ions = QStringLiteral("topNMs2Ions");
    const QString kMinFoundMzPeaks = QStringLiteral("minFoundMzPeaks");
    const QString kDeisotopeScans = QStringLiteral("deisotopeScans");

}


//TODO validate these better.
Err PythiaParameterReader::loadPythiaParameters(PythiaParameters *pythiaParameters) {

    using namespace PythiaParameterReaderConstants;

     const QStringList expectedRaggeds = {
             kNTermRagged,
             kCTermRagged,
             kNoRagged,
             kBothRagged
     };

    ERR_INIT
    e = ErrorUtils::isFalse(jsonContentsIsEmpty()); ree
    e = validateJsonKeys(); ree

    const QMap<QString, QVariant> json = jsonContents();
    for(auto it = json.begin(); it != json.end(); it++){

        const QString &jsonKey = it.key();
        const QVariant &jsonValue = it.value();

        if (jsonKey == kNTermCleavePoints){
            pythiaParameters->nTermCleavePoints = jsonValue.toStringList();
        }
        else if (jsonKey == kCTermCleavePoints){
            pythiaParameters->cTermCleavePoints = jsonValue.toStringList();
        }
        else if(jsonKey == kRaggedness){

            const QString raggedness = jsonValue.toString();
            e = ErrorUtils::isTrue(expectedRaggeds.contains(raggedness)); ree;

            if (raggedness == kNoRagged){
                pythiaParameters->raggedness = Raggedness::NoRagged;
            }
            else if (raggedness == kCTermRagged){
                pythiaParameters->raggedness = Raggedness::CTermRagged;
            }
            else if (raggedness == kNTermRagged){
                pythiaParameters->raggedness = Raggedness::NTermRagged;
            }
            else if (raggedness == kBothRagged){
                pythiaParameters->raggedness = Raggedness::BothRagged;
            }

        }
        else if (jsonKey == kAllowedMissedCleavages){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->allowedMissedCleavages = val;
        }
        else if (jsonKey == kMzMinDataStructure){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->mzMinDataStructure = val;
        }
        else if (jsonKey == kMzMaxDataStructure){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->mzMaxDataStructure = val;
        }
        else if (jsonKey == kPeptideLengthMin){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->peptideLengthMin = val;
        }
        else if (jsonKey == kPeptideLengthMax){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->peptideLengthMax = val;
        }
        else if (jsonKey == kMaxModificationsPeptide){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->maxModificationsPeptide = val;
        }
        else if (jsonKey == kModifications) {

            const QVariantList modsList = jsonValue.toList();
            for (const QVariant &v : modsList) {

                const QMap<QString, QVariant> m =  v.toMap();

                Modification mod;
                for (const QString &key : Modification::modKeys()) {

                    if (key == Modification::modKeys().at(0)) {

                        const QString modResidueOrPositionalLocation = m.value(key).toString().remove(" ");

                        if (modResidueOrPositionalLocation.size() == 1) {
                            mod.residue = modResidueOrPositionalLocation[0];
                        }
                        else {
                            mod.positionalLocation = modResidueOrPositionalLocation;
                        }
                    }
                    if (key == Modification::modKeys().at(1)) {
                        mod.name = m.value(key).toString();
                    }
                    if (key == Modification::modKeys().at(2)) {
                        mod.type = m.value(key).toString() == kFixed
                                   ? ModificationType::FIXED
                                   : ModificationType::DYNAMIC;
                    }
                    if (key == Modification::modKeys().at(3)) {
                        mod.formula = m.value(key).toString();
                    }
                }

                pythiaParameters->modifications.push_back(mod);
            }
        }
        else if (jsonKey == kPrecursorExtractionWindowThomsons){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->precursorExtractionWindowThomsons = val;
        }
        else if (jsonKey == kPercentFDR){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->percentFDR = val;
        }
        else if (jsonKey == kMaxTandemPointCount){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->maxTandemPointCount = val;
        }
        else if (jsonKey == kReturnPSMTopN){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->returnPSMTopN = val;
        }
        else if (jsonKey == kMS2ExtractionWidthPPM){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->ms2ExtractionWidthPPM = val;
        }
        else if (jsonKey == kChargeStateMin){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->chargeStateMin = val;
        }
        else if (jsonKey == kChargeStateMax){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->chargeStateMax = val;
        }
        else if (jsonKey == kAddDecoys){
            bool val = jsonValue.toBool();
            pythiaParameters->addDecoys = val;
        }
        else if (jsonKey == kDeisotopeScans){
            bool val = jsonValue.toBool();
            pythiaParameters->deisotopeScans = val;
        }
        else if (jsonKey == kTopNMs2Ions){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->topNMs2Ions = val;
        }
        else if (jsonKey == kSkipScanCount){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->skipScanCount = val;
        }
        else if (jsonKey == kMinScanCount){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->minScanCount = val;
        }
        else if (jsonKey == kUseMeanMz){
            bool val = jsonValue.toBool();
            pythiaParameters->useMeanMz = val;
        }
        else if (jsonKey == kFilterLength){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->filterLength = val;
        }
        else if (jsonKey == kSmoothCount){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->smoothCount = val;
        }
        else if (jsonKey == kSigma){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->sigma = val;
        }
        else if (jsonKey == kSignalToNoiseRatio){
            double val;
            e = ErrorUtils::toDouble(jsonValue, &val); ree;
            pythiaParameters->signalToNoiseRatio = val;
        }
        else if (jsonKey == kMinFoundMzPeaks){
            int val;
            e = ErrorUtils::toInt(jsonValue, &val); ree;
            pythiaParameters->minFoundMzPeaks = val;
        }
    }

    e = applyFixedModificationsToAminoAcids(
            *pythiaParameters,
            &pythiaParameters->aminoAcids
            );ree

    ERR_RETURN
}


Err PythiaParameterReader::validateJsonKeys() {

    using namespace PythiaParameterReaderConstants;

    ERR_INIT

    const QStringList expectedJsonKeys = {
            kNTermCleavePoints,
            kCTermCleavePoints,
            kRaggedness,
            kAllowedMissedCleavages,
            kMzMinDataStructure,
            kMzMaxDataStructure,
            kPeptideLengthMin,
            kPeptideLengthMax,
            kModifications,
            kMaxModificationsPeptide,

            kMaxTandemPointCount,
            kReturnPSMTopN,
            kChargeStateMin,
            kChargeStateMax,
            kMS2ExtractionWidthPPM,
            kPrecursorExtractionWindowThomsons,
            kPercentFDR,
            kAddDecoys,

            kTopNMs2Ions,
            kSkipScanCount,
            kMinScanCount,
            kUseMeanMz,
            kFilterLength,
            kSmoothCount,
            kSigma,
            kSignalToNoiseRatio,
            kMinFoundMzPeaks
    };

    //TODO FIX THIS
//    const QStringList json = jsonContents().keys();
//    e = ErrorUtils::isEqual(json.size(), expectedJsonKeys.size()); ree;
//
//    for(const QString &key : json){
//        e = ErrorUtils::isTrue(json.contains(key)); ree;
//    }

    ERR_RETURN
}

Err PythiaParameterReader::applyFixedModificationsToAminoAcids(
        const PythiaParameters &parameters,
        AminoAcids *aminoAcids
        ) {

    ERR_INIT

    for (const Modification &mod : parameters.modifications) {

        if (mod.type == ModificationType::DYNAMIC) {
            continue;
        }

        MolecularFormula mf;
        e = parseMolecularFormulaString(mod.formula, &mf); ree;
        aminoAcids->addFixedModification(mod.residue, mf);
    }

    ERR_RETURN
}

PythiaParameters PythiaParameterReader::genericPythiaParametersForTests() {

    PythiaParameters pythiaParameters;

    pythiaParameters.returnPSMTopN = 500;
    pythiaParameters.maxTandemPointCount = 2;
    pythiaParameters.ms2ExtractionWidthPPM = 17.0;
    pythiaParameters.precursorExtractionWindowThomsons = 0.0;
    pythiaParameters.chargeStateMin = 2;
    pythiaParameters.chargeStateMax = 3;
    pythiaParameters.minScanCount = 2;
    pythiaParameters.skipScanCount = 0;
    pythiaParameters.useMeanMz = true;
    pythiaParameters.filterLength = 5;
    pythiaParameters.smoothCount = 2;
    pythiaParameters.sigma = 1;
    pythiaParameters.signalToNoiseRatio = 2;
    pythiaParameters.topNMs2Ions = 6;
    pythiaParameters.minFoundMzPeaks = 3;
    pythiaParameters.allowedMissedCleavages = 1;
    pythiaParameters.mzMinDataStructure = 176.0;
    pythiaParameters.mzMaxDataStructure = 1500.0;
    pythiaParameters.pValThreshold = 0.05;

    Modification carboxyAmidoMethyl(
            'C',
            "CAM",
            ModificationType::FIXED,
            "C2H3NO"
    );
    pythiaParameters.modifications = {carboxyAmidoMethyl};

    PythiaParameterReader::applyFixedModificationsToAminoAcids(
            pythiaParameters,
            &pythiaParameters.aminoAcids
    );

    return pythiaParameters;
}

Err PythiaParameterReader::buildPythiaParameters(
        const QString &pythiaParametersFilePath,
        PythiaParameters *pythiaParameters
) {

    ERR_INIT

    PythiaParameterReader pythiaParameterReader;
    e = pythiaParameterReader.readFile(pythiaParametersFilePath); ree;

    e = pythiaParameterReader.loadPythiaParameters(pythiaParameters); ree;
    pythiaParameters->print();

    ERR_RETURN
}
