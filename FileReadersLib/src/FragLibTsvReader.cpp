//
// Created by anichols on 2/20/24.
//

#include "FragLibTsvReader.h"

#include "ChemConstants.h"
#include "FragLibReader.h"
#include "PeptideStringWithMods.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


namespace {

    std::vector<std::string> splitOnTab(
            const std::string &line,
            int expectedSize
            ) {

        std::vector<std::string> result(expectedSize);
        std::istringstream ss(line);
        std::string field;

        int i = 0;
        while (std::getline(ss, field, '\t')) {
            result[i++] = field;
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

//    Err addIonLabelLogic(FragLibReaderRow *fragLibReaderRow) {
//
//        ERR_INIT
//
//        PeptideSequenceChargeKey peptideSequenceChargeKey = fragLibReaderRow->peptideSequenceChargeKey;
//
//        const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
//                S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
//                Qt::SkipEmptyParts
//                );
//
//        e = ErrorUtils::isEqual(peptideSequenceChargeKeySplit.size(), 2); ree;
//
//        int charge;
//        e = ErrorUtils::toInt(peptideSequenceChargeKeySplit.back(), &charge); ree;
//        const PeptideStringWithMods peptideStringWithMods = PeptideStringWithMods(peptideSequenceChargeKeySplit.front());
//
//
//
//        ERR_RETURN
//    }
//
//    Err addIonLabelInformation(QVector<FragLibReaderRow> *fragLibReaderRows) {
//
//        ERR_INIT
//
//        for (FragLibReaderRow &flrr : *fragLibReaderRows) {
//            e = addIonLabelLogic(&flrr); ree;
//        }
//
//
//        ERR_RETURN
//    }

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
    QList<int> headerIndexes;

    std::ifstream file(fragLibFilePath.toStdString());
    std::string line;

    int columnCount = -1;

    if (file) {
        while(std::getline(file, line)){

            if (columnCount < 0) {
                columnCount = static_cast<int>(std::count(line.begin(), line.end(), '\t')) + 1;
            }

            const std::vector<std::string> lineSplit = splitOnTab(line, columnCount);

            if (header.empty()) {
                header = lineSplit;
                headerIndexVsColumnNames = buildHeaderIndexVsColumnNames(lineSplit);
                headerIndexes = headerIndexVsColumnNames.keys();
                e = ErrorUtils::isNotEmpty(headerIndexVsColumnNames); ree;
                continue;
            }

            FragLibTsvReaderRow fragLibTsvReaderRow;

            for (int headerIndex : headerIndexes) {

                const std::string &colName = headerIndexVsColumnNames.value(headerIndex);
                if (colName.empty()) {
                    continue;
                }

                const std::string &valString = lineSplit.at(headerIndex);
                if (valString.empty()) {
                    continue;
                }

                if (colName == PRECURSOR_MZ) {
                    fragLibTsvReaderRow.precursorMz = std::stod(valString);
//                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.precursorMz); ree;
                }
                else if (colName == PRODUCT_MZ){
                    fragLibTsvReaderRow.productMz = std::stod(valString);
//                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.productMz); ree;
                }
                else if (colName == TR_RECALIB){
                    fragLibTsvReaderRow.trRecalibrated = std::stod(valString);
//                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.trRecalibrated); ree;
                }
                else if (colName == ION_MOBILITY){
                    fragLibTsvReaderRow.ionMobility = std::stod(valString);
//                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.ionMobility); ree;
                }
                else if (colName == MOD_PEP){
                    fragLibTsvReaderRow.modifiedPeptide = QString::fromStdString(valString);
                }
                else if (colName == PRECURSOR_CHARGE){
                    fragLibTsvReaderRow.precursorCharge = std::stoi(valString);
//                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.precursorCharge); ree;
                }
                else if (colName == LIB_INTENSITY){
                    fragLibTsvReaderRow.libraryIntensity = std::stod(valString);
//                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.libraryIntensity); ree;
                }
                else if (colName == DECOY) {
                    fragLibTsvReaderRow.decoy = std::stoi(valString);
//                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.decoy); ree;
                }
                else if (colName == FRAG_TYPE){
                    fragLibTsvReaderRow.fragmentType = QString::fromStdString(valString);
                }
                else if (colName == FRAG_CHARGE){
                    fragLibTsvReaderRow.fragmentCharge = std::stoi(valString);
//                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentCharge); ree;
                }
                else if (colName == FRAG_SERIES_NUMB){
                    fragLibTsvReaderRow.fragmentSeriesNumber = std::stoi(valString);
//                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentSeriesNumber); ree;
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

    m_fragLibReaderRows.swap(*fragLibReaderRows);

    if (!m_mzValsCurrent.isEmpty()) {
        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.precursorCharge = m_precursorChargeCurrent;
        fragLibReaderRow.iRT = m_irtCurrent;
        fragLibReaderRow.intensityVals = m_intensityCurrent;
        fragLibReaderRow.mzVals = m_mzValsCurrent;
        fragLibReaderRow.ionLabels = m_labelsCurrent.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.isDecoy = m_decoyCurrent;
        fragLibReaderRow.peptideSequenceChargeKey = m_currentPeptide;
        fragLibReaderRow.mass = (m_precursorMzCurrent * static_cast<float>(m_precursorChargeCurrent))
                                - (ChemConstants::PROTON * m_precursorChargeCurrent);

        if (massStart <= fragLibReaderRow.mass && fragLibReaderRow.mass <= massEnd) {
            fragLibReaderRows->push_back(fragLibReaderRow);
        }

    }

//    if (fragLibReaderRows->front().ionLabels.contains("-1^")) {
//        addIonLabelInformation(fragLibReaderRows);
//    }

    qDebug() << "MS2 Predictions count:" << fragLibReaderRows->size() << "retrieved in" << et.elapsed() << "mSec";
    ERR_RETURN


}

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

        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.precursorCharge = m_precursorChargeCurrent;
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
    if (tsvRow.fragmentCharge > 0) {
        m_labelsCurrent.push_back(ionLabel);
    }
    m_mzValsCurrent.push_back(static_cast<float>(tsvRow.productMz));
    m_intensityCurrent.push_back(static_cast<float>(tsvRow.libraryIntensity));
    m_irtCurrent = static_cast<float>(tsvRow.trRecalibrated);
    m_precursorChargeCurrent = tsvRow.precursorCharge;
    m_precursorMzCurrent = static_cast<float>(tsvRow.precursorMz);
    m_decoyCurrent = tsvRow.decoy;

    ERR_RETURN
}
