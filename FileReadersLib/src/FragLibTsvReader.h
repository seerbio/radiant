//
// Created by anichols on 2/20/24.
//

#ifndef PYTHIADIACPP_FRAGLIBTSVREADER_H
#define PYTHIADIACPP_FRAGLIBTSVREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "PeptideStringWithMods.h"

#include <QHash>

using namespace Error;


namespace FragLibTsvReaderRowNamespace {
    /* DIANN tsv columns:
        FileName
        PrecursorMz
        ProductMz
        Tr_recalibrated
        IonMobility
        transition_name
        LibraryIntensity
        transition_group_id
        decoy
        PeptideSequence
        Proteotypic
        QValue
        PGQValue
        Ms1ProfileCorr
        ProteinGroup
        ProteinName
        Genes
        FullUniModPeptideName
        ModifiedPeptide
        PrecursorCharge
        PeptideGroupLabel
        UniprotID
        NTerm
        CTerm
        FragmentType
        FragmentCharge
        FragmentSeriesNumber
        FragmentLossType
        ExcludeFromAssay
    */

    enum class FragLibTsvColumn {
        PrecursorMz,
        ProductMz,
        TrRecalibrated,
        IonMobility,
        LibraryIntensity,
        Decoy,
        ModifiedPeptide,
        PrecursorCharge,
        FragmentType,
        FragmentCharge,
        FragmentSeriesNumber
    };

    struct FragLibTsvColumnBinding {
        int headerIndex = -1;
        FragLibTsvColumn column = FragLibTsvColumn::PrecursorMz;
    };

    struct FragLibTsvColumnAliasPriority {
        FragLibTsvColumn column = FragLibTsvColumn::PrecursorMz;
        QVector<QString> aliases;
    };

    const QString FRAG_LIB_TSV_SUFFIX = QStringLiteral("tsv");

    const QString PRECURSOR_MZ = QStringLiteral("PrecursorMz");
    const QString PRODUCT_MZ = QStringLiteral("ProductMz");
    const QString FRAGMENT_MZ = QStringLiteral("FragmentMz");
    const QString TR_RECALIB = QStringLiteral("Tr_recalibrated");
    const QString RT = QStringLiteral("RT");
    const QString NORM_RT = QStringLiteral("NormalizedRetentionTime");
    const QString ION_MOBILITY = QStringLiteral("IonMobility");
    const QString LIB_INTENSITY = QStringLiteral("LibraryIntensity");
    const QString RELATIVE_INTENSITY = QStringLiteral("RelativeIntensity");
    const QString DECOY = QStringLiteral("decoy");
    const QString MOD_PEP = QStringLiteral("ModifiedPeptide");
    const QString MOD_PEP_SEQ = QStringLiteral("ModifiedPeptideSequence");
    const QString PRECURSOR_CHARGE = QStringLiteral("PrecursorCharge");
    const QString FRAG_TYPE = QStringLiteral("FragmentType");
    const QString FRAG_CHARGE = QStringLiteral("FragmentCharge");
    const QString FRAG_SERIES_NUMB = QStringLiteral("FragmentSeriesNumber");
    const QString FRAG_NUMB = QStringLiteral("FragmentNumber");

    inline const QVector<FragLibTsvColumnAliasPriority> &columnAliases() {
        static const QVector<FragLibTsvColumnAliasPriority> stacks = {
                {FragLibTsvColumn::PrecursorMz, {PRECURSOR_MZ}},
                {FragLibTsvColumn::ProductMz, {PRODUCT_MZ, FRAGMENT_MZ}},
                {FragLibTsvColumn::TrRecalibrated, {TR_RECALIB, RT, NORM_RT}},
                {FragLibTsvColumn::IonMobility, {ION_MOBILITY}},
                {FragLibTsvColumn::LibraryIntensity, {LIB_INTENSITY, RELATIVE_INTENSITY}},
                {FragLibTsvColumn::Decoy, {DECOY}},
                {FragLibTsvColumn::ModifiedPeptide, {MOD_PEP, MOD_PEP_SEQ}},
                {FragLibTsvColumn::PrecursorCharge, {PRECURSOR_CHARGE}},
                {FragLibTsvColumn::FragmentType, {FRAG_TYPE}},
                {FragLibTsvColumn::FragmentCharge, {FRAG_CHARGE}},
                {FragLibTsvColumn::FragmentSeriesNumber, {FRAG_SERIES_NUMB, FRAG_NUMB}}
        };

        return stacks;
    }

}//namespace

/**
 * Purposely not using CSVReaderInputBase because it is much slower for processing millions of lines. due to QMap building
 */
struct FragLibTsvReaderRow {

    double precursorMz = -1.0;
    double productMz = -1.0;
    double trRecalibrated = -1.0;
    double ionMobility = -1.0;
    double libraryIntensity = -1.0;
    int decoy = -1;
    QString modifiedPeptide;
    int precursorCharge = -1;
    QString fragmentType;
    int fragmentCharge = -1;
    int fragmentSeriesNumber = -1;

};


class FragLibReaderRow;

class FILEREADERSLIB_EXPORTS FragLibTsvReader {

public:

    // Err writeToFragLibFormat();

    Err getFragLibReaderRows(
            const QString &fragLibFilePath,
            QList<FragLibReaderRow> *fragLibReaderRows
    );

    void setEnableTerminalByPenultimateDecoyAnnotationShift(bool enable);

    static Err inferIonLabelsForTest(
            const QVector<float> &mzValsToPair,
            const PeptideStringWithMods &peptideStringWithMods,
            bool isDecoy,
            int charge,
            bool enableTerminalByPenultimateDecoyAnnotationShift,
            QString *ionLabels
            );

private:

    Err convertFragLibTsvReaderRowsToFragLibReaderRow(const FragLibTsvReaderRow &tsvRow);

private:

    QList<FragLibReaderRow> m_fragLibReaderRows;
    QString m_tsvFilePath;

    QStringList m_labelsCurrent;
    QVector<float> m_mzValsCurrent;
    QVector<float> m_intensityCurrent;
    QString m_currentPeptide;
    int m_decoyCurrent = -1;
    float m_irtCurrent = 1.0;
    int m_precursorChargeCurrent = -1;
    float m_precursorMzCurrent = -1.0;
    bool m_enableTerminalByPenultimateDecoyAnnotationShift = false;

};


#endif //PYTHIADIACPP_FRAGLIBTSVREADER_H
