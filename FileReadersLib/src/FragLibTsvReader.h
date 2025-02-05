//
// Created by anichols on 2/20/24.
//

#ifndef PYTHIADIACPP_FRAGLIBTSVREADER_H
#define PYTHIADIACPP_FRAGLIBTSVREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"


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

    const QString FRAG_LIB_TSV_SUFFIX = QStringLiteral("tsv");

    const QString PRECURSOR_MZ = "PrecursorMz";
    const QString PRODUCT_MZ = "ProductMz";
    const QString TR_RECALIB = "Tr_recalibrated";
    const QString ION_MOBILITY = "IonMobility";
    const QString LIB_INTENSITY = "LibraryIntensity";
    const QString DECOY = "decoy";
    const QString MOD_PEP = "ModifiedPeptide";
    const QString PRECURSOR_CHARGE = "PrecursorCharge";
    const QString FRAG_TYPE = "FragmentType";
    const QString FRAG_CHARGE = "FragmentCharge";
    const QString FRAG_SERIES_NUMB = "FragmentSeriesNumber";

    const QVector<QString> keysToCheck = {
            PRECURSOR_MZ,
            PRODUCT_MZ,
            TR_RECALIB,
            ION_MOBILITY,
            MOD_PEP,
            PRECURSOR_CHARGE,
            LIB_INTENSITY,
            DECOY,
            FRAG_TYPE,
            FRAG_CHARGE,
            FRAG_SERIES_NUMB

    };
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
            QVector<FragLibReaderRow> *fragLibReaderRows
    );

private:

    Err convertFragLibTsvReaderRowsToFragLibReaderRow(const FragLibTsvReaderRow &tsvRow);

private:

    QVector<FragLibReaderRow> m_fragLibReaderRows;
    QString m_tsvFilePath;

    QStringList m_labelsCurrent;
    QVector<float> m_mzValsCurrent;
    QVector<float> m_intensityCurrent;
    QString m_currentPeptide;
    int m_decoyCurrent = -1;
    float m_irtCurrent = 1.0;
    int m_precursorChargeCurrent = -1;
    float m_precursorMzCurrent = -1.0;

};


#endif //PYTHIADIACPP_FRAGLIBTSVREADER_H
