//
// Created by Drucifer on 12/4/2021.
//

#include "PythiaParameterReader.h"
#include "ErrorUtils.h"

#include "toml.hpp"

#include <QDebug>
#include <QJsonValue>

#include <iostream>
#include <sstream>


Modification::Modification(
        QChar residue,
        const QString &name,
        const ModificationType &type,
        const QString &formula
        )
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
    const QString kNTermCleavePoints = QStringLiteral("nTermCleavePoints");
    const QString kCTermCleavePoints = QStringLiteral("cTermCleavePoints");
    const QString kRaggedness = QStringLiteral("raggedness");
    const QString kAllowedMissedCleavages = QStringLiteral("allowedMissedCleavages");
    const QString kFixed = QStringLiteral("fixed");
    const QString kDynamic = QStringLiteral("dynamic");
    const QString kMaxModificationsPeptide = QStringLiteral("maxModificationsPeptide");
    const QString kModifications = QStringLiteral("modifications");
    const QString kMzMinMS2  = QStringLiteral("mzMinMS2");
    const QString kMzMaxMS2  = QStringLiteral("mzMaxMS2");
    const QString kPeptideLengthMin = QStringLiteral("peptideLengthMin");
    const QString kPeptideLengthMax = QStringLiteral("peptideLengthMax");
    const QString kNoRagged = QStringLiteral("NoRagged");
    const QString kCTermRagged = QStringLiteral("CTermRagged");
    const QString kNTermRagged = QStringLiteral("NTermRagged");
    const QString kBothRagged = QStringLiteral("BothRagged");

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
    const QString kMinFoundMzPeaks = QStringLiteral("minFoundMzPeaks");
    const QString kCosineSimToAnchorThreshold = QStringLiteral("cosineSimToAnchorThreshold");
    const QString kScanTimeWindowStDevs = QStringLiteral("scanTimeWindowStDevs");
    const QString kReportDecoys = QStringLiteral("reportDecoys");
    const QString  kSubtractShadows = QStringLiteral("subtractShadows");

    const QString kDigestParams = QStringLiteral("DigestParams");
    const QString kMS2Params = QStringLiteral("MS2Params");
    const QString kPrecursorParams = QStringLiteral("PrecursorParams");
    const QString kPeakIntegrationParams = QStringLiteral("PeakIntegrationParams");
    const QString kFdrParams = QStringLiteral("FdrParams");
}

Err PythiaParameterReader::validateJsonKeys() {

    using namespace PythiaParameterReaderConstants;

    ERR_INIT

    const QStringList expectedJsonKeys = {
            kNTermCleavePoints,
            kCTermCleavePoints,
            kRaggedness,
            kAllowedMissedCleavages,
            kMzMinMS2,
            kMzMaxMS2,
            kPeptideLengthMin,
            kPeptideLengthMax,
            kModifications,
            kMaxModificationsPeptide,

            kChargeStateMin,
            kChargeStateMax,
            kMS2ExtractionWidthPPM,
            kPrecursorExtractionWindowThomsons,
            kPercentFDR,

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
    pythiaParameters.minFoundMzPeaks = 3;
    pythiaParameters.allowedMissedCleavages = 1;
    pythiaParameters.mzMinMS2 = 176.0;
    pythiaParameters.mzMaxMS2 = 1500.0;

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

    //TODO validate these better.

    using namespace PythiaParameterReaderConstants;

    e = ErrorUtils::fileExists(pythiaParametersFilePath); ree;
    auto parser = toml::parse_file(pythiaParametersFilePath.toStdString());

    const auto digestNode =  parser[kDigestParams.toStdString()];

    toml::array* cTermCleavePoints = digestNode[kCTermCleavePoints.toStdString()].as_array();
    for(const auto& value : *cTermCleavePoints) {
        pythiaParameters->cTermCleavePoints.push_back(value.value_or(""));
    }

    toml::array* nTermCleavePoints = digestNode[kNTermCleavePoints.toStdString()].as_array();
    for(const auto& value : *nTermCleavePoints) {
        pythiaParameters->nTermCleavePoints.push_back(value.value_or(""));
    }

    const QString raggedness = QString::fromStdString(std::string(digestNode[kRaggedness.toStdString()].value_or("")));
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

    pythiaParameters->allowedMissedCleavages = digestNode[kAllowedMissedCleavages.toStdString()].value_or(0);
    pythiaParameters->peptideLengthMin = digestNode[kPeptideLengthMin.toStdString()].value_or(0);
    pythiaParameters->peptideLengthMax = digestNode[kPeptideLengthMax.toStdString()].value_or(0);
    pythiaParameters->maxModificationsPeptide = digestNode[kMaxModificationsPeptide.toStdString()].value_or(0);

    const auto ms2ParamsNode =  parser[kMS2Params.toStdString()];
    pythiaParameters->minFoundMzPeaks = ms2ParamsNode[kMinFoundMzPeaks.toStdString()].value_or(0);
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

    toml::array* modifications = parser["Modification"].as_array();
    for(const auto& value : *modifications) {
        const toml::table *table = value.as_table();

        Modification mod;

        const QString modResidueOrPositionalLocation = QString::fromStdString(std::string(table->at("residue").value_or("")));

        if (modResidueOrPositionalLocation.size() == 1) {
            mod.residue = modResidueOrPositionalLocation[0];
        }
        else {
            mod.positionalLocation = modResidueOrPositionalLocation;
        }

        mod.name = QString::fromStdString(std::string(table->at("name").value_or("")));

        const QString modType = QString::fromStdString(std::string(table->at("type").value_or("")));
        mod.type = modType == kFixed ? ModificationType::FIXED : ModificationType::DYNAMIC;

        mod.formula = QString::fromStdString(std::string(table->at("formula").value_or("")));

        pythiaParameters->modifications.push_back(mod);
    }

    e = applyFixedModificationsToAminoAcids(
            *pythiaParameters,
            &pythiaParameters->aminoAcids
    );ree


    ERR_RETURN
}
