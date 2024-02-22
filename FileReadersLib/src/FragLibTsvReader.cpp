//
// Created by anichols on 2/20/24.
//

#include "FragLibTsvReader.h"

#include "ChemConstants.h"
#include "FragLibReader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


QMap<QString, double> FragLibTsvReader::uniModMap() {

    const QMap<QString, double> uniModNameVsModificationMass = {
            {QStringLiteral("UniMod:4"), static_cast<double>(57.021464)},
            {QStringLiteral("Carbamidomethyl (C)"), static_cast<double>(57.021464)},
            {QStringLiteral("Carbamidomethyl"), static_cast<double>(57.021464)},
            {QStringLiteral("CAM"), static_cast<double>(57.021464)},
            {QStringLiteral("+57"), static_cast<double>(57.021464)},
            {QStringLiteral("+57.0"), static_cast<double>(57.021464)},
            {QStringLiteral("UniMod:26"), static_cast<double>(39.994915)},
            {QStringLiteral("PCm"), static_cast<double>(39.994915)},
            {QStringLiteral("UniMod:5"), static_cast<double>(43.005814)},
            {QStringLiteral("Carbamylation (KR)"), static_cast<double>(43.005814)},
            {QStringLiteral("+43"), static_cast<double>(43.005814)},
            {QStringLiteral("+43.0"), static_cast<double>(43.005814)},
            {QStringLiteral("CRM"), static_cast<double>(43.005814)},
            {QStringLiteral("UniMod:7"), static_cast<double>(0.984016)},
            {QStringLiteral("Deamidation (NQ)"), static_cast<double>(0.984016)},
            {QStringLiteral("Deamidation"), static_cast<double>(0.984016)},
            {QStringLiteral("Dea"), static_cast<double>(0.984016)},
            {QStringLiteral("+1"), static_cast<double>(0.984016)},
            {QStringLiteral("+1.0"), static_cast<double>(0.984016)},
            {QStringLiteral("UniMod:35"), static_cast<double>(15.994915)},
            {QStringLiteral("Oxidation (M)"), static_cast<double>(15.994915)},
            {QStringLiteral("Oxidation"), static_cast<double>(15.994915)},
            {QStringLiteral("Oxi"), static_cast<double>(15.994915)},
            {QStringLiteral("+16"), static_cast<double>(15.994915)},
            {QStringLiteral("+16.0"), static_cast<double>(15.994915)},
            {QStringLiteral("Oxi"), static_cast<double>(15.994915)},
            {QStringLiteral("UniMod:1"), static_cast<double>(42.010565)},
            {QStringLiteral("Acetyl (Protein N-term)"), static_cast<double>(42.010565)},
            {QStringLiteral("+42"), static_cast<double>(42.010565)},
            {QStringLiteral("+42.0"), static_cast<double>(42.010565)},
            {QStringLiteral("UniMod:255"), static_cast<double>(28.0313)},
            {QStringLiteral("AAR"), static_cast<double>(28.0313)},
            {QStringLiteral("UniMod:254"), static_cast<double>(26.01565)},
            {QStringLiteral("AAS"), static_cast<double>(26.01565)},
            {QStringLiteral("UniMod:122"), static_cast<double>(27.994915)},
            {QStringLiteral("Frm"), static_cast<double>(27.994915)},
            {QStringLiteral("UniMod:1301"), static_cast<double>(128.094963)},
            {QStringLiteral("+1K"), static_cast<double>(128.094963)},
            {QStringLiteral("UniMod:1288"), static_cast<double>(156.101111)},
            {QStringLiteral("+1R"), static_cast<double>(156.101111)},
            {QStringLiteral("UniMod:27"), static_cast<double>(-18.010565)},
            {QStringLiteral("PGE"), static_cast<double>(-18.010565)},
            {QStringLiteral("UniMod:28"), static_cast<double>(-17.026549)},
            {QStringLiteral("PGQ"), static_cast<double>(-17.026549)},
            {QStringLiteral("UniMod:526"), static_cast<double>(-48.003371)},
            {QStringLiteral("DTM"), static_cast<double>(-48.003371)},
            {QStringLiteral("UniMod:325"), static_cast<double>(31.989829)},
            {QStringLiteral("2Ox"), static_cast<double>(31.989829)},
            {QStringLiteral("UniMod:342"), static_cast<double>(15.010899)},
            {QStringLiteral("Amn"), static_cast<double>(15.010899)},
            {QStringLiteral("UniMod:1290"), static_cast<double>(114.042927)},
            {QStringLiteral("2CM"), static_cast<double>(114.042927)},
            {QStringLiteral("UniMod:359"), static_cast<double>(13.979265)},
            {QStringLiteral("PGP"), static_cast<double>(13.979265)},
            {QStringLiteral("UniMod:30"), static_cast<double>(21.981943)},
            {QStringLiteral("NaX"), static_cast<double>(21.981943)},
            {QStringLiteral("UniMod:401"), static_cast<double>(-2.015650)},
            {QStringLiteral("-2H"), static_cast<double>(-2.015650)},
            {QStringLiteral("UniMod:528"), static_cast<double>(14.999666)},
            {QStringLiteral("MDe"), static_cast<double>(14.999666)},
            {QStringLiteral("UniMod:385"), static_cast<double>(-17.026549)},
            {QStringLiteral("dAm"), static_cast<double>(-17.026549)},
            {QStringLiteral("UniMod:23"), static_cast<double>(-18.010565)},
            {QStringLiteral("Dhy"), static_cast<double>(-18.010565)},
            {QStringLiteral("UniMod:129"), static_cast<double>(125.896648)},
            {QStringLiteral("Iod"), static_cast<double>(125.896648)},
            {QStringLiteral("Phosphorylation (ST)"), static_cast<double>(79.966331)},
            {QStringLiteral("UniMod:21"), static_cast<double>(79.966331)},
            {QStringLiteral("+80"), static_cast<double>(79.966331)},
            {QStringLiteral("+80.0"), static_cast<double>(79.966331)},
            {QStringLiteral("UniMod:259"), static_cast<double>(8.014199)},
            {QStringLiteral("Lys8"),static_cast<double>(8.014199)},
            {QStringLiteral("UniMod:267"), static_cast<double>(10.008269)},
            {QStringLiteral("Arg10"), static_cast<double>(10.008269)},
            {QStringLiteral("UniMod:268"), static_cast<double>(6.013809)},
            {QStringLiteral("UniMod:269"), static_cast<double>(10.027228)}
    };

    return uniModNameVsModificationMass;
}

namespace {

    std::vector<std::string> splitOnTab(const std::string &line) {

        std::vector<std::string> result;
        std::istringstream ss(line);
        std::string field;

        while (std::getline(ss, field, '\t')) {
            result.push_back(field);
        }

        return result;
    }

    QHash<int, std::string> buildHeaderIndexVsColumnNames(const std::vector<std::string> &headerLine) {

        QHash<int, std::string> hashTable;

        int counter = 0;
        for (const std::string &ss : headerLine) {

            if (FragLibTsvReaderRowNamespace::keysToCheck.contains(ss)) {
                hashTable.insert(counter, ss);
            }

            counter++;
        }

        return hashTable;
    }

}//namespace
Err FragLibTsvReader::getFragLibReaderRows(
        const QString &fragLibFilePath,
        double massStart,
        double massEnd,
        QVector<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    using namespace FragLibTsvReaderRowNamespace;

    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(massEnd > massStart); ree;

    m_fragLibReaderRows.clear();

    QFileInfo fi(fragLibFilePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            FragLibTsvReaderRowNamespace::FRAG_LIB_TSV_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

    fragLibReaderRows->clear();

    QElapsedTimer et;
    et.start();

    std::vector<std::string> header;
    QHash<int, std::string> headerIndexVsColumnNames;
    const QList<int> &headerIndexes = headerIndexVsColumnNames.keys();

    std::ifstream file(fragLibFilePath.toStdString());
    std::string line;

    if (file) {
        while(std::getline(file, line)){

            const std::vector<std::string> lineSplit = splitOnTab(line);

            if (header.empty()) {
                header = lineSplit;
                headerIndexVsColumnNames = buildHeaderIndexVsColumnNames(lineSplit);
                e = ErrorUtils::isNotEmpty(headerIndexVsColumnNames); ree;
                continue;
            }

            FragLibTsvReaderRow fragLibTsvReaderRow;

            for (int headerIndex : headerIndexes) {

                const std::string &colName = headerIndexVsColumnNames.value(headerIndex);
                if (colName.empty()) {
                    continue;
                }

                const QString valString = QString::fromStdString(lineSplit.at(headerIndex));

                if (colName == PRECURSOR_MZ) {
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.precursorMz); ree;
                }
                else if (colName == PRODUCT_MZ){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.productMz); ree;
                }
                else if (colName == TR_RECALIB){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.trRecalibrated); ree;
                }
                else if (colName == ION_MOBILITY){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.ionMobility); ree;
                }
                else if (colName == MOD_PEP){
                    fragLibTsvReaderRow.modifiedPeptide = valString;
                }
                else if (colName == PRECURSOR_CHARGE){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.precursorCharge); ree;
                }
                else if (colName == LIB_INTENSITY){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.libraryIntensity); ree;
                }
                else if (colName == DECOY) {
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.decoy); ree;
                }
                else if (colName == FRAG_TYPE){
                    fragLibTsvReaderRow.fragmentType = valString;
                }
                else if (colName == FRAG_CHARGE){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentCharge); ree;
                }
                else if (colName == FRAG_SERIES_NUMB){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentSeriesNumber); ree;
                }
                else {
                    rrr(eValueError);
                }
            }

            e = convertFragLibTsvReaderRowsToFragLibReaderRow(
                    fragLibTsvReaderRow,
                    massStart,
                    massEnd
            ); ree;
        }
    }

    if (!m_mzValsCurrent.isEmpty()) {
        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.iRT = m_irtCurrent;
        fragLibReaderRow.intensityVals = m_intensityCurrent;
        fragLibReaderRow.mzVals = m_mzValsCurrent;
        fragLibReaderRow.ionLabels = m_labelsCurrent.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.isDecoy = m_decoyCurrent;
        fragLibReaderRow.peptideSequenceChargeKey = m_currentPeptide;
        fragLibReaderRow.mass = (m_precursorMzCurrent * static_cast<float>(m_precursorChargeCurrent))
                                - (ChemConstants::PROTON * m_precursorChargeCurrent);

        if (massStart <= fragLibReaderRow.mass && fragLibReaderRow.mass <= massEnd) {
            m_fragLibReaderRows.push_back(fragLibReaderRow);
        }
    }

    qDebug() << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN


}

static int counterr = 0;
Err FragLibTsvReader::convertFragLibTsvReaderRowsToFragLibReaderRow(
        const FragLibTsvReaderRow &tsvRow,
        double massStart,
        double massEnd
) {

    ERR_INIT

    const QString peptideSequenceChargeKey = tsvRow.modifiedPeptide + "|" + QString::number(tsvRow.precursorCharge);
    if (m_currentPeptide.isEmpty()) {
        m_currentPeptide = peptideSequenceChargeKey;
    }

    if (m_currentPeptide != peptideSequenceChargeKey) {

        std::cout << ++counterr << std::endl;

        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.iRT = m_irtCurrent;
        fragLibReaderRow.intensityVals = m_intensityCurrent;
        fragLibReaderRow.mzVals = m_mzValsCurrent;
        fragLibReaderRow.ionLabels = m_labelsCurrent.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.isDecoy = m_decoyCurrent;
        fragLibReaderRow.peptideSequenceChargeKey = m_currentPeptide;
        fragLibReaderRow.mass = (m_precursorMzCurrent * static_cast<float>(m_precursorChargeCurrent))
                                - (ChemConstants::PROTON * m_precursorChargeCurrent);

        if (massStart <= fragLibReaderRow.mass && fragLibReaderRow.mass <= massEnd) {
            m_fragLibReaderRows.push_back(fragLibReaderRow);
        }

        m_intensityCurrent.clear();
        m_mzValsCurrent.clear();
        m_labelsCurrent.clear();
        m_currentPeptide = peptideSequenceChargeKey;
    }

    QString ionLabel = tsvRow.fragmentType + QString::number(tsvRow.fragmentSeriesNumber);
    ionLabel += tsvRow.fragmentCharge == 1 ? "" : "^" + QString::number(tsvRow.fragmentCharge);
    m_labelsCurrent.push_back(ionLabel);
    m_mzValsCurrent.push_back(static_cast<float>(tsvRow.productMz));
    m_intensityCurrent.push_back(static_cast<float>(tsvRow.libraryIntensity));
    m_irtCurrent = static_cast<float>(tsvRow.trRecalibrated);
    m_precursorChargeCurrent = tsvRow.precursorCharge;
    m_precursorMzCurrent = static_cast<float>(tsvRow.precursorMz);

    ERR_RETURN
}
