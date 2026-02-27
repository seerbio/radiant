//
// Created by anichols on 2/20/24.
//

#include "FragLibTsvReader.h"

#include "AminoAcids.h"
#include "ChemConstants.h"
#include "FragLibReader.h"
#include "FragLibReaderRow.h"
#include "PeptideStringWithMods.h"

#include <nanoflann.hpp>
#include <Eigen/Dense>

#include <QtConcurrent/QtConcurrent>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


namespace {
    void maybeShiftSeriesForDecoyAnnotation(
            const bool isDecoy,
            const bool enableTerminalByPenultimateDecoyAnnotationShift,
            const bool isYSeries,
            const int ionCharge,
            const double nTermDeltaMass,
            const double cTermDeltaMass,
            QVector<double> *seriesMzVals
            ) {

        if (!isDecoy || !enableTerminalByPenultimateDecoyAnnotationShift) {
            return;
        }

        const double chargeScale = ionCharge > 1 ? (1.0 / static_cast<double>(ionCharge)) : 1.0;
        for (int i = 0; i < seriesMzVals->size(); ++i) {
            const int ionIndex = i + 1;
            const int residueIndexCoveredByFragment = ionIndex; // terminal-by-penultimate mode

            double mzShift = 0.0;
            if (isYSeries) {
                if (residueIndexCoveredByFragment >= 1) {
                    mzShift += cTermDeltaMass * chargeScale;
                }
                if (residueIndexCoveredByFragment >= seriesMzVals->size() - 1) {
                    mzShift += nTermDeltaMass * chargeScale;
                }
            }
            else {
                if (residueIndexCoveredByFragment >= 1) {
                    mzShift += nTermDeltaMass * chargeScale;
                }
                if (residueIndexCoveredByFragment >= seriesMzVals->size() - 1) {
                    mzShift += cTermDeltaMass * chargeScale;
                }
            }

            (*seriesMzVals)[i] += mzShift;
        }
    }

    QVector<QString> splitOnTab(
            const std::string &line,
            int expectedSize
            ) {

        QVector<QString> result;
        result.reserve(expectedSize);

        QString lineQ = QString::fromStdString(line);
        lineQ = lineQ.replace("\r", "").replace("\r\n", "");

        QStringList fields = lineQ.split('\t');

        for(const QString &field : fields) {
            result.push_back(field);
        }

        return result;
    }

    QHash<int, QString> buildHeaderIndexVsColumnNames(const QVector<QString> &headerLine) {

        QHash<int, QString> hashTable;

        int counter = 0;
        for (const QString &ss : headerLine) {
            if (FragLibTsvReaderRowNamespace::keysToCheck.contains(ss)) {
                hashTable.insert(counter, ss);
            }

            counter++;
        }

        return hashTable;
    }

    Err splitPeptideSequenceChargeKey(
            const FragLibReaderRow &fragLibReaderRow,
            PeptideStringWithMods *peptideStringWithMods,
            int *charge
    ) {

        ERR_INIT

        PeptideSequenceChargeKey peptideSequenceChargeKey = fragLibReaderRow.peptideSequenceChargeKey;

        const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                Qt::SkipEmptyParts
        );

        e = ErrorUtils::isEqual(peptideSequenceChargeKeySplit.size(), 2); ree;
        e = ErrorUtils::toInt(peptideSequenceChargeKeySplit.back(), charge); ree;
        *peptideStringWithMods = PeptideStringWithMods(peptideSequenceChargeKeySplit.front());

        ERR_RETURN
    }

    using KDTree = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd>;

    Err loadKDTreeStructures(
        const PeptideStringWithMods &peptideStringWithMods,
        bool isDecoy,
        bool enableTerminalByPenultimateDecoyAnnotationShift,
        int charge,
        QVector<double> *mzValsTree,
        QStringList *ionLabelsList,
        Eigen::MatrixX<double> *mat
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;

        mzValsTree->clear();
        ionLabelsList->clear();

        AminoAcids aminoAcids;
        double nTermDeltaMass = 0.0;
        double cTermDeltaMass = 0.0;
        const PeptideString peptideString = peptideStringWithMods.removeUniModChars();
        if (peptideString.size() >= 3) {
            constexpr int firstIndexToMutate = 1;
            const int secondIndexToMutate = peptideString.size() - 2;
            const QMap<QChar, double> diannMutateAminoAcidToMass = AminoAcids::diannMutateAminoAcidToMass();
            nTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(firstIndexToMutate));
            cTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(secondIndexToMutate));
        }

        QVector<double> bSeriesMz = peptideStringWithMods.bSeries(1, aminoAcids);
        maybeShiftSeriesForDecoyAnnotation(
                isDecoy,
                enableTerminalByPenultimateDecoyAnnotationShift,
                false,
                1,
                nTermDeltaMass,
                cTermDeltaMass,
                &bSeriesMz
                );
        mzValsTree->append(bSeriesMz);
        ionLabelsList->append(peptideStringWithMods.bSeriesIonLabels({}));

        QVector<double> ySeriesMz = peptideStringWithMods.ySeries(1, aminoAcids);
        maybeShiftSeriesForDecoyAnnotation(
                isDecoy,
                enableTerminalByPenultimateDecoyAnnotationShift,
                true,
                1,
                nTermDeltaMass,
                cTermDeltaMass,
                &ySeriesMz
                );
        mzValsTree->append(ySeriesMz);
        ionLabelsList->append(peptideStringWithMods.ySeriesIonLabels({}));

        if (charge > 1) {

            for (int i = 1; i < charge; i++) {

                const QString chargeChar = QString::number(i+1);

                QVector<double> bSeriesMzN = peptideStringWithMods.bSeries(i + 1, aminoAcids);
                maybeShiftSeriesForDecoyAnnotation(
                        isDecoy,
                        enableTerminalByPenultimateDecoyAnnotationShift,
                        false,
                        i + 1,
                        nTermDeltaMass,
                        cTermDeltaMass,
                        &bSeriesMzN
                        );
                mzValsTree->append(bSeriesMzN);
                ionLabelsList->append(peptideStringWithMods.bSeriesIonLabels("^" + chargeChar));

                QVector<double> ySeriesMzN = peptideStringWithMods.ySeries(i + 1, aminoAcids);
                maybeShiftSeriesForDecoyAnnotation(
                        isDecoy,
                        enableTerminalByPenultimateDecoyAnnotationShift,
                        true,
                        i + 1,
                        nTermDeltaMass,
                        cTermDeltaMass,
                        &ySeriesMzN
                        );
                mzValsTree->append(ySeriesMzN);
                ionLabelsList->append(peptideStringWithMods.ySeriesIonLabels("^" + chargeChar));
            }
        }

        e = ErrorUtils::isEqual(mzValsTree->size(), ionLabelsList->size());

        mat->resize(mzValsTree->size() ,2);
        mat->setZero();
        for (int i = 0; i < mzValsTree->size(); i++) {
            mat->coeffRef(i, 0) = mzValsTree->at(i);
            mat->coeffRef(i, 1) = 0.0;
        }

        ERR_RETURN
    }

    Err matchMzValsToIonType(
            const QVector<float> &mzValsToPair,
            const PeptideStringWithMods &peptideStringWithMods,
            bool isDecoy,
            int charge,
            bool enableTerminalByPenultimateDecoyAnnotationShift,
            QString *ionLabels
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzValsToPair); ree;
        e = ErrorUtils::isNotEmpty(peptideStringWithMods); ree;

        QVector<double> mzValsTree;
        QStringList ionLabelsList;
        Eigen::MatrixX<double> mat;
        e = loadKDTreeStructures(
            peptideStringWithMods,
            isDecoy,
            enableTerminalByPenultimateDecoyAnnotationShift,
            charge,
            &mzValsTree,
            &ionLabelsList,
            &mat
            ); ree;

        const int treeLeafSize = 5;
        KDTree kdTree(static_cast<int>(mat.cols()), mat, treeLeafSize);

        QStringList ionLabelsFoundList;
        for (float mzVal : mzValsToPair) {

            if (mzVal < 1) {
                continue;
            }

            const size_t numResults = 1;
            std::vector<double> queryPt = {mzVal, 0.0};
            std::vector<long> mzIndex(numResults);
            std::vector<double> outDistSqr(numResults);

            std::vector<std::pair<Eigen::Index, double>> matches;

            const size_t resultsSize = kdTree.index->knnSearch(
                    queryPt.data(),
                    numResults,
                    mzIndex.data(),
                    outDistSqr.data()
            );
            e = ErrorUtils::isTrue(resultsSize > 0); ree;

            const double mzValClosestFound = mzValsTree.at(static_cast<int>(mzIndex.front()));
            const double distanceFromFound = std::abs(mzValClosestFound - mzVal);

            if (distanceFromFound > 0.2) {
                qDebug() << peptideStringWithMods;
                qDebug() << "Decoy:" << isDecoy;
                qDebug() << "Searched mzVal:" << mzVal;
                qDebug() << "Closest val found:" << mzValClosestFound;
                qDebug() << "Distance from closest val found:" << distanceFromFound;
                qDebug() << "Theo Frag" << mzValsTree;
                qDebug() << "Emp Frag" << mzValsToPair;
                rrr(eValueError);
            }

            ionLabelsFoundList.push_back(ionLabelsList.value(static_cast<int>(mzIndex.front())));
        }

        *ionLabels = ionLabelsFoundList.join(S_GLOBAL_SETTINGS.SEPARATOR);

        ERR_RETURN
    }

    Err buildIonLabels(
            const FragLibReaderRow &fragLibReaderRow,
            bool enableTerminalByPenultimateDecoyAnnotationShift,
            QString *ionLabels
    ) {

        ERR_INIT

        PeptideStringWithMods peptideStringWithMods;
        int charge;
        e = splitPeptideSequenceChargeKey(
                fragLibReaderRow,
                &peptideStringWithMods,
                &charge
        ); ree;

        e = matchMzValsToIonType(
                fragLibReaderRow.mzVals,
                peptideStringWithMods,
                fragLibReaderRow.isDecoy,
                charge,
                enableTerminalByPenultimateDecoyAnnotationShift,
                ionLabels
        ); ree;

        ERR_RETURN
    }

    void parallelLogic(
            FragLibReaderRow &fragLibReaderRow,
            bool enableTerminalByPenultimateDecoyAnnotationShift
            ) {

        ERR_INIT

        e = buildIonLabels(
                fragLibReaderRow,
                enableTerminalByPenultimateDecoyAnnotationShift,
                &fragLibReaderRow.ionLabels
                );

        if (e != eNoError) {
            throw std::runtime_error("Building Ion Labels didn't work out"); einfo;
        }

    }

    Err addIonLabelInformation(
            QList<FragLibReaderRow> *fragLibReaderRows,
            bool enableTerminalByPenultimateDecoyAnnotationShift
            ) {

        ERR_INIT

#define PARALLEL_LABELS
#ifdef PARALLEL_LABELS
    const auto buildLabelsLogic = [enableTerminalByPenultimateDecoyAnnotationShift](FragLibReaderRow &fragLibReaderRow) {
        parallelLogic(fragLibReaderRow, enableTerminalByPenultimateDecoyAnnotationShift);
    };
    QFuture<void> futures = QtConcurrent::map(*fragLibReaderRows, buildLabelsLogic);
    futures.waitForFinished();
#else
        for (FragLibReaderRow &flrr : *fragLibReaderRows) {
            parallelLogic(flrr, enableTerminalByPenultimateDecoyAnnotationShift);
        }
#endif

        ERR_RETURN
    }

    Err normalizeIntensityValues(
        const QVector<float> &intensityVals,
        QVector<float> *intensityValsNormalized
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(intensityVals); ree;
        intensityValsNormalized->clear();

        const float intensityMax = *std::max_element(intensityVals.begin(), intensityVals.end());
        std::transform(
            intensityVals.begin(),
            intensityVals.end(),
            std::back_inserter(*intensityValsNormalized),
            [intensityMax](const float val) {return val / std::max(1.0f, intensityMax); }
            );

        ERR_RETURN
    }

}//namespace
Err FragLibTsvReader::getFragLibReaderRows(
        const QString &tsvFilePath,
        QList<FragLibReaderRow> *fragLibReaderRows
        ) {

    ERR_INIT

    using namespace FragLibTsvReaderRowNamespace;

    e = ErrorUtils::fileExists(tsvFilePath); ree;

    m_tsvFilePath = tsvFilePath;

    m_fragLibReaderRows.clear();

    QFileInfo fi(tsvFilePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            FragLibTsvReaderRowNamespace::FRAG_LIB_TSV_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

    fragLibReaderRows->clear();

    QElapsedTimer et;
    et.start();

    QVector<QString> header;
    QHash<int, QString> headerIndexVsColumnNames;
    QList<int> headerIndexes;

    std::ifstream file(tsvFilePath.toStdString());
    std::string line;

    int columnCount = -1;

    if (file) {
        while(std::getline(file, line)){

            if (columnCount < 0) {
                columnCount = static_cast<int>(std::count(line.begin(), line.end(), '\t')) + 1;
            }

            const QVector<QString> lineSplit = splitOnTab(line, columnCount);

            if (header.empty()) {
                header = lineSplit;
                headerIndexVsColumnNames = buildHeaderIndexVsColumnNames(lineSplit);
                headerIndexes = headerIndexVsColumnNames.keys();
                e = ErrorUtils::isNotEmpty(headerIndexVsColumnNames); ree;
                continue;
            }

            FragLibTsvReaderRow fragLibTsvReaderRow;

            for (int headerIndex : headerIndexes) {

                const QString &colName = headerIndexVsColumnNames.value(headerIndex);

                if (colName.isEmpty()) {
                    continue;
                }

                const QString &valString = lineSplit.at(headerIndex);
                if (valString.isEmpty()) {
                    continue;
                }

                if (colName == PRECURSOR_MZ) {
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.precursorMz); eee_absorb;
                }
                else if (colName == PRODUCT_MZ){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.productMz); eee_absorb;
                }
                else if (colName == TR_RECALIB){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.trRecalibrated); eee_absorb;
                }
                else if (colName == ION_MOBILITY){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.ionMobility); eee_absorb;
                }
                else if (colName == MOD_PEP || colName == MOD_PEP_SEQ){
                    fragLibTsvReaderRow.modifiedPeptide = valString;
                }
                else if (colName == PRECURSOR_CHARGE){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.precursorCharge); eee_absorb;
                }
                else if (colName == LIB_INTENSITY){
                    e = ErrorUtils::toDouble(valString, &fragLibTsvReaderRow.libraryIntensity); eee_absorb;
                }
                else if (colName == DECOY) {
                	const QString boolConversion = valString.toLower() == "true" || valString == "1" ? "1" : "0";
                    e = ErrorUtils::toInt(boolConversion, &fragLibTsvReaderRow.decoy); eee_absorb;
                }
                else if (colName == FRAG_TYPE){
                    fragLibTsvReaderRow.fragmentType = valString;
                }
                else if (colName == FRAG_CHARGE){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentCharge); eee_absorb;
                }
                else if (colName == FRAG_SERIES_NUMB){
                    e = ErrorUtils::toInt(valString, &fragLibTsvReaderRow.fragmentSeriesNumber); eee_absorb;
                }
                else {
                    rrr(eValueError);
                }
            }

            e = convertFragLibTsvReaderRowsToFragLibReaderRow(fragLibTsvReaderRow); ree;
        }
    }

    m_fragLibReaderRows.swap(*fragLibReaderRows);

    if (!m_mzValsCurrent.isEmpty()) {

        QVector<float> intensityCurrentNormalized;
        e = normalizeIntensityValues(
            m_intensityCurrent,
            &intensityCurrentNormalized
            ); ree;

        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.precursorCharge = m_precursorChargeCurrent;
        fragLibReaderRow.iRT = m_irtCurrent;
        fragLibReaderRow.intensityVals = intensityCurrentNormalized;
        fragLibReaderRow.mzVals = m_mzValsCurrent;
        fragLibReaderRow.ionLabels = m_labelsCurrent.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.isDecoy = m_decoyCurrent;
        fragLibReaderRow.peptideSequenceChargeKey = m_currentPeptide;
        fragLibReaderRow.mass = (m_precursorMzCurrent * static_cast<float>(m_precursorChargeCurrent))
                                - (ChemConstants::PROTON * m_precursorChargeCurrent);

        fragLibReaderRows->push_back(fragLibReaderRow);
    }

    if (fragLibReaderRows->front().ionLabels.contains("-1^") || fragLibReaderRows->front().ionLabels.isEmpty()) {
        e = addIonLabelInformation(
                fragLibReaderRows,
                m_enableTerminalByPenultimateDecoyAnnotationShift
                ); ree;
    }

    ERR_RETURN
}

void FragLibTsvReader::setEnableTerminalByPenultimateDecoyAnnotationShift(bool enable) {
    m_enableTerminalByPenultimateDecoyAnnotationShift = enable;
}

Err FragLibTsvReader::inferIonLabelsForTest(
        const QVector<float> &mzValsToPair,
        const PeptideStringWithMods &peptideStringWithMods,
        bool isDecoy,
        int charge,
        bool enableTerminalByPenultimateDecoyAnnotationShift,
        QString *ionLabels
        ) {

    return matchMzValsToIonType(
            mzValsToPair,
            peptideStringWithMods,
            isDecoy,
            charge,
            enableTerminalByPenultimateDecoyAnnotationShift,
            ionLabels
            );
}

Err FragLibTsvReader::convertFragLibTsvReaderRowsToFragLibReaderRow(
    const FragLibTsvReaderRow &tsvRow
    ) {

    ERR_INIT
    const QString peptideSequenceChargeKey = tsvRow.modifiedPeptide + "|" + QString::number(tsvRow.precursorCharge);
    if (m_currentPeptide.isEmpty()) {
        m_currentPeptide = peptideSequenceChargeKey;
    }

    if (m_currentPeptide != peptideSequenceChargeKey) {

        QVector<float> intensityValsNormalized;
        e = normalizeIntensityValues(
                   m_intensityCurrent,
                   &intensityValsNormalized
                   ); ree;

        FragLibReaderRow fragLibReaderRow;
        fragLibReaderRow.precursorCharge = m_precursorChargeCurrent;
        fragLibReaderRow.iRT = m_irtCurrent;
        fragLibReaderRow.intensityVals = intensityValsNormalized;
        fragLibReaderRow.mzVals = m_mzValsCurrent;
        fragLibReaderRow.ionLabels = m_labelsCurrent.join(S_GLOBAL_SETTINGS.SEPARATOR);
        fragLibReaderRow.isDecoy = m_decoyCurrent;
        fragLibReaderRow.peptideSequenceChargeKey = m_currentPeptide;
        fragLibReaderRow.mass = (m_precursorMzCurrent * static_cast<float>(m_precursorChargeCurrent))
                                - (ChemConstants::PROTON * m_precursorChargeCurrent);

        m_fragLibReaderRows.push_back(fragLibReaderRow);

        m_intensityCurrent.clear();
        m_mzValsCurrent.clear();
        m_labelsCurrent.clear();
        m_currentPeptide = peptideSequenceChargeKey;
    }

    if (tsvRow.productMz < 1.0) {
        ERR_RETURN
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
