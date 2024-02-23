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
