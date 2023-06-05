//
// Created by anichols on 5/30/23.
//

#ifndef PYTHIADIACPP_FRAMEEXTRACTSREADER_H
#define PYTHIADIACPP_FRAMEEXTRACTSREADER_H

#include <utility>

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

using namespace Error;


namespace FrameExtractsReaderNamespace {

    const QString PEP_WITH_MODS= QStringLiteral("peptideStringWithMods");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString SCAN_NUMBER = QStringLiteral("scanNumber");
    const QString FRAME_INDEX = QStringLiteral("frameIndex");
    const QString COSINE_SUM_WEIGHTED = QStringLiteral("cosineSimSumWeighted");

    const QString A_IND_THEO = QStringLiteral("aIonIndexesTheo");
    const QString A_MZ_THEO = QStringLiteral("aIonMzValsTheo");
    const QString A_INTZ_THEO = QStringLiteral("aIonIntesitiesTheo");
    const QString A_IND_ACT = QStringLiteral("aIonIndexesActual");
    const QString A_MZ_ACT = QStringLiteral("aIonMzValsActual");
    const QString A_INTZ_ACT = QStringLiteral("aIonIntesitiesActual");
    const QString A_HILL_COS_SIM = QStringLiteral("aIonHillCosineSims");
    const QString Y_IND_THEO = QStringLiteral("yIonIndexesTheo");
    const QString Y_MZ_THEO = QStringLiteral("yIonMzValsTheo");
    const QString Y_INTZ_THEO = QStringLiteral("yIonIntesitiesTheo");
    const QString Y_IND_ACT = QStringLiteral("yIonIndexesActual");
    const QString Y_MZ_ACT = QStringLiteral("yIonMzValsActual");
    const QString Y_INTZ_ACT = QStringLiteral("yIonIntesitiesActual");
    const QString Y_HILL_COS_SIM = QStringLiteral("yIonHillCosineSims");
    const QString Y2_IND_THEO = QStringLiteral("y2IonIndexesTheo");
    const QString Y2_MZ_THEO = QStringLiteral("y2IonMzValsTheo");
    const QString Y2_INTZ_THEO = QStringLiteral("y2IonIntesitiesTheo");
    const QString Y2_IND_ACT = QStringLiteral("y2IonIndexesActual");
    const QString Y2_MZ_ACT = QStringLiteral("y2IonMzValsActual");
    const QString Y2_INTZ_ACT = QStringLiteral("y2IonIntesitiesActual");
    const QString Y2_HILL_COS_SIM = QStringLiteral("y2IonHillCosineSims");
    const QString YNH3_IND_THEO = QStringLiteral("yNH3IonIndexesTheo");
    const QString YNH3_MZ_THEO = QStringLiteral("yNH3IonMzValsTheo");
    const QString YNH3_INTZ_THEO = QStringLiteral("yNH3IonIntesitiesTheo");
    const QString YNH3_IND_ACT = QStringLiteral("yNH3IonIndexesActual");
    const QString YNH3_MZ_ACT = QStringLiteral("yNH3IonMzValsActual");
    const QString YNH3_INTZ_ACT = QStringLiteral("yNH3IonIntesitiesActual");
    const QString YNH3_HILL_COS_SIM = QStringLiteral("yNH3IonHillCosineSims");
    const QString YH2O_IND_THEO = QStringLiteral("yH2OIonIndexesTheo");
    const QString YH2O_MZ_THEO = QStringLiteral("yH2OIonMzValsTheo");
    const QString YH2O_INTZ_THEO = QStringLiteral("yH2OIonIntesitiesTheo");
    const QString YH2O_IND_ACT = QStringLiteral("yH2OIonIndexesActual");
    const QString YH2O_MZ_ACT = QStringLiteral("yH2OIonMzValsActual");
    const QString YH2O_INTZ_ACT = QStringLiteral("yH2OIonIntesitiesActual");
    const QString YH2O_HILL_COS_SIM = QStringLiteral("yH2OIonHillCosineSims");
    const QString B_IND_THEO = QStringLiteral("bIonIndexesTheo");
    const QString B_MZ_THEO = QStringLiteral("bIonMzValsTheo");
    const QString B_INTZ_THEO = QStringLiteral("bIonIntesitiesTheo");
    const QString B_IND_ACT = QStringLiteral("bIonIndexesActual");
    const QString B_MZ_ACT = QStringLiteral("bIonMzValsActual");
    const QString B_INTZ_ACT = QStringLiteral("bIonIntesitiesActual");
    const QString B_HILL_COS_SIM = QStringLiteral("bIonHillCosineSims");
    const QString B2_IND_THEO = QStringLiteral("b2IonIndexesTheo");
    const QString B2_MZ_THEO = QStringLiteral("b2IonMzValsTheo");
    const QString B2_INTZ_THEO = QStringLiteral("b2IonIntesitiesTheo");
    const QString B2_IND_ACT = QStringLiteral("b2IonIndexesActual");
    const QString B2_MZ_ACT = QStringLiteral("b2IonMzValsActual");
    const QString B2_INTZ_ACT = QStringLiteral("b2IonIntesitiesActual");
    const QString B2_HILL_COS_SIM = QStringLiteral("b2IonHillCosineSims");
    const QString BNH3_IND_THEO = QStringLiteral("bNH3IonIndexesTheo");
    const QString BNH3_MZ_THEO = QStringLiteral("bNH3IonMzValsTheo");
    const QString BNH3_INTZ_THEO = QStringLiteral("bNH3IonIntesitiesTheo");
    const QString BNH3_IND_ACT = QStringLiteral("bNH3IonIndexesActual");
    const QString BNH3_MZ_ACT = QStringLiteral("bNH3IonMzValsActual");
    const QString BNH3_INTZ_ACT = QStringLiteral("bNH3IonIntesitiesActual");
    const QString BNH3_HILL_COS_SIM = QStringLiteral("bNH3IonHillCosineSims");
    const QString BH2O_IND_THEO = QStringLiteral("bH2OIonIndexesTheo");
    const QString BH2O_MZ_THEO = QStringLiteral("bH2OIonMzValsTheo");
    const QString BH2O_INTZ_THEO = QStringLiteral("bH2OIonIntesitiesTheo");
    const QString BH2O_IND_ACT = QStringLiteral("bH2OIonIndexesActual");
    const QString BH2O_MZ_ACT = QStringLiteral("bH2OIonMzValsActual");
    const QString BH2O_INTZ_ACT = QStringLiteral("bH2OIonIntesitiesActual");
    const QString BH2O_HILL_COS_SIM = QStringLiteral("bH2OIonHillCosineSims");
    const QString PREC_IND_THEO = QStringLiteral("precursorIonIndexesTheo");
    const QString PREC_MZ_THEO= QStringLiteral("precursorIonMzValsTheo");
    const QString PREC_INTZ_THEO = QStringLiteral("precursorIonIntesitiesTheo");
    const QString PREC_IND_ACT = QStringLiteral("precursorIonIndexesActual");
    const QString PREC_MZ_ACT = QStringLiteral("precursorIonMzValsActual");
    const QString PREC_INTZ_ACT = QStringLiteral("precursorIonIntesitiesActual");
    const QString PREC_HILL_COS_SIM = QStringLiteral("precursorIonHillCosineSims");

    const QString Y_MZ_STD_ACT = QStringLiteral("yIonMzStdActual");
    const QString Y_HILL_LEN_ACT = QStringLiteral("yIonHillLengthActual");
    const QString Y2_MZ_STD_ACT = QStringLiteral("y2IonMzStdActual");
    const QString Y2_HILL_LEN_ACT = QStringLiteral("y2IonHillLengthActual");
    const QString YNH3_MZ_STD_ACT = QStringLiteral("yNH3IonMzStdActual");
    const QString YNH3_HILL_LEN_ACT = QStringLiteral("yNH3IonHillLengthActual");
    const QString YH2O_MZ_STD_ACT = QStringLiteral("yH2OIonMzStdActual");
    const QString YH2O_HILL_LEN_ACT = QStringLiteral("yH2OIonHillLengthActual");
    const QString B_MZ_STD_ACT = QStringLiteral("bIonMzStdActual");
    const QString B_HILL_LEN_ACT = QStringLiteral("bIonHillLengthActual");
    const QString B2_MZ_STD_ACT = QStringLiteral("b2IonMzStdActual");
    const QString B2_HILL_LEN_ACT = QStringLiteral("b2IonHillLengthActual");
    const QString BNH3_MZ_STD_ACT = QStringLiteral("bNH3IonMzStdActual");
    const QString BNH3_HILL_LEN_ACT = QStringLiteral("bNH3IonHillLengthActual");
    const QString BH2O_MZ_STD_ACT = QStringLiteral("bH2OIonMzStdActual");
    const QString BH2O_HILL_LEN_ACT = QStringLiteral("bH2OIonHillLengthActual");
    const QString A_MZ_STD_ACT = QStringLiteral("aIonMzStdActual");
    const QString A_HILL_LEN_ACT = QStringLiteral("aIonHillLengthActual");
    const QString PREC_MZ_STD_ACT = QStringLiteral("precursorIonMzStdActual");
    const QString PREC_HILL_LEN_ACT = QStringLiteral("precursorIonHillLengthActual");

    const QString A_ISO_COS_ACT = QStringLiteral("aIonIsotopologueCosineSimActual");
    const QString A_ISO_CHRG_ACT = QStringLiteral("aIonIsotopologueChargeActual");
    const QString A_ISO_INTZ_ACT = QStringLiteral("aIonIsotopologueIntensityActual");
    const QString Y_ISO_COS_ACT = QStringLiteral("yIonIsotopologueCosineSimActual");
    const QString Y_ISO_CHRG_ACT = QStringLiteral("yIonIsotopologueChargeActual");
    const QString Y_ISO_INTZ_ACT = QStringLiteral("yIonIsotopologueIntensityActual");
    const QString Y2_ISO_COS_ACT = QStringLiteral("y2IonIsotopologueCosineSimActual");
    const QString Y2_ISO_CHRG_ACT = QStringLiteral("y2IonIsotopologueChargeActual");
    const QString Y2_ISO_INTZ_ACT = QStringLiteral("y2IonIsotopologueIntensityActual");
    const QString YNH3_ISO_COS_ACT = QStringLiteral("yNH3IonIsotopologueCosineSimActual");
    const QString YNH3_ISO_CHRG_ACT = QStringLiteral("yNH3IonIsotopologueChargeActual");
    const QString YNH3_ISO_INTZ_ACT = QStringLiteral("yNH3IonIsotopologueIntensityActual");
    const QString YH2O_ISO_COS_ACT = QStringLiteral("yH2OIonIsotopologueCosineSimActual");
    const QString YH2O_ISO_CHRG_ACT = QStringLiteral("yH2OIonIsotopologueChargeActual");
    const QString YH2O_ISO_INTZ_ACT = QStringLiteral("yH2OIonIsotopologueIntensityActual");
    const QString B_ISO_COS_ACT = QStringLiteral("bIonIsotopologueCosineSimActual");
    const QString B_ISO_CHRG_ACT = QStringLiteral("bIonIsotopologueChargeActual");
    const QString B_ISO_INTZ_ACT = QStringLiteral("bIonIsotopologueIntensityActual");
    const QString B2_ISO_COS_ACT = QStringLiteral("b2IonIsotopologueCosineSimActual");
    const QString B2_ISO_CHRG_ACT = QStringLiteral("b2IonIsotopologueChargeActual");
    const QString B2_ISO_INTZ_ACT = QStringLiteral("b2IonIsotopologueIntensityActual");
    const QString BNH3_ISO_COS_ACT = QStringLiteral("bNH3IonIsotopologueCosineSimActual");
    const QString BNH3_ISO_CHRG_ACT = QStringLiteral("bNH3IonIsotopologueChargeActual");
    const QString BNH3_ISO_INTZ_ACT = QStringLiteral("bNH3IonIsotopologueIntensityActual");
    const QString BH2O_ISO_COS_ACT = QStringLiteral("bH2OIonIsotopologueCosineSimActual");
    const QString BH2O_ISO_CHRG_ACT = QStringLiteral("bH2OIonIsotopologueChargeActual");
    const QString BH2O_ISO_INTZ_ACT = QStringLiteral("bH2OIonIsotopologueIntensityActual");
    const QString PREC_ISO_COS_ACT = QStringLiteral("precursorIonIsotopologueCosineSimActual");
    const QString PREC_ISO_CHRG_ACT = QStringLiteral("precursorIonIsotopologueChargeActual");
    const QString PREC_ISO_INTZ_ACT = QStringLiteral("precursorIonIsotopologueIntensityActual");

    const QStringList keysToCheck = {
            PEP_WITH_MODS,
            IS_DECOY,
            SCAN_NUMBER,
            FRAME_INDEX,
            COSINE_SUM_WEIGHTED,
            A_IND_THEO,
            A_MZ_THEO,
            A_INTZ_THEO,
            A_IND_ACT,
            A_MZ_ACT,
            A_INTZ_ACT,
            A_HILL_COS_SIM,
            Y_IND_THEO,
            Y_MZ_THEO,
            Y_INTZ_THEO,
            Y_IND_ACT,
            Y_MZ_ACT,
            Y_INTZ_ACT,
            Y_HILL_COS_SIM,
            Y2_IND_THEO,
            Y2_MZ_THEO,
            Y2_INTZ_THEO,
            Y2_IND_ACT,
            Y2_MZ_ACT,
            Y2_INTZ_ACT,
            Y2_HILL_COS_SIM,
            YNH3_IND_THEO,
            YNH3_MZ_THEO,
            YNH3_INTZ_THEO,
            YNH3_IND_ACT,
            YNH3_MZ_ACT,
            YNH3_INTZ_ACT,
            YNH3_HILL_COS_SIM,
            YH2O_IND_THEO,
            YH2O_MZ_THEO,
            YH2O_INTZ_THEO,
            YH2O_IND_ACT,
            YH2O_MZ_ACT,
            YH2O_INTZ_ACT,
            YH2O_HILL_COS_SIM,
            B_IND_THEO,
            B_MZ_THEO,
            B_INTZ_THEO,
            B_IND_ACT,
            B_MZ_ACT,
            B_INTZ_ACT,
            B_HILL_COS_SIM,
            B2_IND_THEO,
            B2_MZ_THEO,
            B2_INTZ_THEO,
            B2_IND_ACT,
            B2_MZ_ACT,
            B2_INTZ_ACT,
            B2_HILL_COS_SIM,
            BNH3_IND_THEO,
            BNH3_MZ_THEO,
            BNH3_INTZ_THEO,
            BNH3_IND_ACT,
            BNH3_MZ_ACT,
            BNH3_INTZ_ACT,
            BNH3_HILL_COS_SIM,
            BH2O_IND_THEO,
            BH2O_MZ_THEO,
            BH2O_INTZ_THEO,
            BH2O_IND_ACT,
            BH2O_MZ_ACT,
            BH2O_INTZ_ACT,
            BH2O_HILL_COS_SIM,
            PREC_IND_THEO,
            PREC_MZ_THEO,
            PREC_INTZ_THEO,
            PREC_IND_ACT,
            PREC_MZ_ACT,
            PREC_INTZ_ACT,
            PREC_HILL_COS_SIM,

            Y_MZ_STD_ACT,
            Y_HILL_LEN_ACT,
            Y2_MZ_STD_ACT,
            Y2_HILL_LEN_ACT,
            YNH3_MZ_STD_ACT,
            YNH3_HILL_LEN_ACT,
            YH2O_MZ_STD_ACT,
            YH2O_HILL_LEN_ACT,
            B_MZ_STD_ACT,
            B_HILL_LEN_ACT,
            B2_MZ_STD_ACT,
            B2_HILL_LEN_ACT,
            BNH3_MZ_STD_ACT,
            BNH3_HILL_LEN_ACT,
            BH2O_MZ_STD_ACT,
            BH2O_HILL_LEN_ACT,
            A_MZ_STD_ACT,
            A_HILL_LEN_ACT,
            PREC_MZ_STD_ACT,
            PREC_HILL_LEN_ACT,

            A_ISO_COS_ACT,
            A_ISO_CHRG_ACT,
            A_ISO_INTZ_ACT,
            Y_ISO_COS_ACT,
            Y_ISO_CHRG_ACT,
            Y_ISO_INTZ_ACT,
            Y2_ISO_COS_ACT,
            Y2_ISO_CHRG_ACT,
            Y2_ISO_INTZ_ACT,
            YNH3_ISO_COS_ACT,
            YNH3_ISO_CHRG_ACT,
            YNH3_ISO_INTZ_ACT,
            YH2O_ISO_COS_ACT,
            YH2O_ISO_CHRG_ACT,
            YH2O_ISO_INTZ_ACT,
            B_ISO_COS_ACT,
            B_ISO_CHRG_ACT,
            B_ISO_INTZ_ACT,
            B2_ISO_COS_ACT,
            B2_ISO_CHRG_ACT,
            B2_ISO_INTZ_ACT,
            BNH3_ISO_COS_ACT,
            BNH3_ISO_CHRG_ACT,
            BNH3_ISO_INTZ_ACT,
            BH2O_ISO_COS_ACT,
            BH2O_ISO_CHRG_ACT,
            BH2O_ISO_INTZ_ACT,
            PREC_ISO_COS_ACT,
            PREC_ISO_CHRG_ACT,
            PREC_ISO_INTZ_ACT
    };

}


struct  FILEREADERSLIB_EXPORTS FrameExtractsReaderRow : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    bool isDecoy = false;
    ScanNumber scanNumberApex = -1;
    FrameIndex frameIndexApex = -1;
    double cosineSimSumWeighted = -1.0;

    QVector<IonIndex> aIonIndexesTheo;
    QVector<double> aIonMzValsTheo;
    QVector<double> aIonIntesitiesTheo;
    QVector<IonIndex> aIonIndexesActual;
    QVector<double> aIonMzValsActual;
    QVector<double> aIonIntesitiesActual;
    QVector<double> aIonHillCosineSims;

    QVector<IonIndex> yIonIndexesTheo;
    QVector<double> yIonMzValsTheo;
    QVector<double> yIonIntesitiesTheo;
    QVector<IonIndex> yIonIndexesActual;
    QVector<double> yIonMzValsActual;
    QVector<double> yIonIntesitiesActual;
    QVector<double> yIonHillCosineSims;

    QVector<IonIndex> y2IonIndexesTheo;
    QVector<double> y2IonMzValsTheo;
    QVector<double> y2IonIntesitiesTheo;
    QVector<IonIndex> y2IonIndexesActual;
    QVector<double> y2IonMzValsActual;
    QVector<double> y2IonIntesitiesActual;
    QVector<double> y2IonHillCosineSims;

    QVector<IonIndex> yNH3IonIndexesTheo;
    QVector<double> yNH3IonMzValsTheo;
    QVector<double> yNH3IonIntesitiesTheo;
    QVector<IonIndex> yNH3IonIndexesActual;
    QVector<double> yNH3IonMzValsActual;
    QVector<double> yNH3IonIntesitiesActual;
    QVector<double> yNH3IonHillCosineSims;

    QVector<IonIndex> yH2OIonIndexesTheo;
    QVector<double> yH2OIonMzValsTheo;
    QVector<double> yH2OIonIntesitiesTheo;
    QVector<IonIndex> yH2OIonIndexesActual;
    QVector<double> yH2OIonMzValsActual;
    QVector<double> yH2OIonIntesitiesActual;
    QVector<double> yH2OIonHillCosineSims;

    QVector<IonIndex> bIonIndexesTheo;
    QVector<double> bIonMzValsTheo;
    QVector<double> bIonIntesitiesTheo;
    QVector<IonIndex> bIonIndexesActual;
    QVector<double> bIonMzValsActual;
    QVector<double> bIonIntesitiesActual;
    QVector<double> bIonHillCosineSims;

    QVector<IonIndex> b2IonIndexesTheo;
    QVector<double> b2IonMzValsTheo;
    QVector<double> b2IonIntesitiesTheo;
    QVector<IonIndex> b2IonIndexesActual;
    QVector<double> b2IonMzValsActual;
    QVector<double> b2IonIntesitiesActual;
    QVector<double> b2IonHillCosineSims;

    QVector<IonIndex> bNH3IonIndexesTheo;
    QVector<double> bNH3IonMzValsTheo;
    QVector<double> bNH3IonIntesitiesTheo;
    QVector<IonIndex> bNH3IonIndexesActual;
    QVector<double> bNH3IonMzValsActual;
    QVector<double> bNH3IonIntesitiesActual;
    QVector<double> bNH3IonHillCosineSims;

    QVector<IonIndex> bH2OIonIndexesTheo;
    QVector<double> bH2OIonMzValsTheo;
    QVector<double> bH2OIonIntesitiesTheo;
    QVector<IonIndex> bH2OIonIndexesActual;
    QVector<double> bH2OIonMzValsActual;
    QVector<double> bH2OIonIntesitiesActual;
    QVector<double> bH2OIonHillCosineSims;

    QVector<IonIndex> precursorIonIndexesTheo;
    QVector<double> precursorIonMzValsTheo;
    QVector<double> precursorIonIntesitiesTheo;
    QVector<IonIndex> precursorIonIndexesActual;
    QVector<double> precursorIonMzValsActual;
    QVector<double> precursorIonIntesitiesActual;
    QVector<double> precursorIonHillCosineSims;

    QVector<double> aIonMzStdActual;
    QVector<int> aIonHillLengthActual;
    QVector<double> yIonMzStdActual;
    QVector<int> yIonHillLengthActual;
    QVector<double> y2IonMzStdActual;
    QVector<int> y2IonHillLengthActual;
    QVector<double> yNH3IonMzStdActual;
    QVector<int> yNH3IonHillLengthActual;
    QVector<double> yH2OIonMzStdActual;
    QVector<int> yH2OIonHillLengthActual;
    QVector<double> bIonMzStdActual;
    QVector<int> bIonHillLengthActual;
    QVector<double> b2IonMzStdActual;
    QVector<int> b2IonHillLengthActual;
    QVector<double> bNH3IonMzStdActual;
    QVector<int> bNH3IonHillLengthActual;
    QVector<double> bH2OIonMzStdActual;
    QVector<int> bH2OIonHillLengthActual;
    QVector<double> precursorIonMzStdActual;
    QVector<int> precursorIonHillLengthActual;

    QVector<double> aIonIsotopologueCosineSimActual;
    QVector<int> aIonIsotopologueChargeActual;
    QVector<double> aIonIsotopologueIntensityActual;
    QVector<double> yIonIsotopologueCosineSimActual;
    QVector<int> yIonIsotopologueChargeActual;
    QVector<double> yIonIsotopologueIntensityActual;
    QVector<double> y2IonIsotopologueCosineSimActual;
    QVector<int> y2IonIsotopologueChargeActual;
    QVector<double> y2IonIsotopologueIntensityActual;
    QVector<double> yNH3IonIsotopologueCosineSimActual;
    QVector<int> yNH3IonIsotopologueChargeActual;
    QVector<double> yNH3IonIsotopologueIntensityActual;
    QVector<double> yH2OIonIsotopologueCosineSimActual;
    QVector<int> yH2OIonIsotopologueChargeActual;
    QVector<double> yH2OIonIsotopologueIntensityActual;
    QVector<double> bIonIsotopologueCosineSimActual;
    QVector<int> bIonIsotopologueChargeActual;
    QVector<double> bIonIsotopologueIntensityActual;
    QVector<double> b2IonIsotopologueCosineSimActual;
    QVector<int> b2IonIsotopologueChargeActual;
    QVector<double> b2IonIsotopologueIntensityActual;
    QVector<double> bNH3IonIsotopologueCosineSimActual;
    QVector<int> bNH3IonIsotopologueChargeActual;
    QVector<double> bNH3IonIsotopologueIntensityActual;
    QVector<double> bH2OIonIsotopologueCosineSimActual;
    QVector<int> bH2OIonIsotopologueChargeActual;
    QVector<double> bH2OIonIsotopologueIntensityActual;
    QVector<double> precursorIonIsotopologueCosineSimActual;
    QVector<int> precursorIonIsotopologueChargeActual;
    QVector<double> precursorIonIsotopologueIntensityActual;


    QMap<QString, QVariant> map() override {

        using namespace FrameExtractsReaderNamespace;

        return {
                {PEP_WITH_MODS, QVariant(peptideStringWithMods)},
                {IS_DECOY, QVariant(isDecoy)},
                {SCAN_NUMBER, QVariant(scanNumberApex)},
                {FRAME_INDEX, QVariant(frameIndexApex)},
                {COSINE_SUM_WEIGHTED, QVariant(cosineSimSumWeighted)},

                {A_IND_THEO, QVariant(qVectorToQByteArray(aIonIndexesTheo))},
                {A_MZ_THEO, QVariant(qVectorToQByteArray(aIonMzValsTheo))},
                {A_INTZ_THEO, QVariant(qVectorToQByteArray(aIonIntesitiesTheo))},
                {A_IND_ACT, QVariant(qVectorToQByteArray(aIonIndexesActual))},
                {A_MZ_ACT, QVariant(qVectorToQByteArray(aIonMzValsActual))},
                {A_INTZ_ACT, QVariant(qVectorToQByteArray(aIonIntesitiesActual))},
                {A_HILL_COS_SIM, QVariant(qVectorToQByteArray(aIonHillCosineSims))},

                {Y_IND_THEO, QVariant(qVectorToQByteArray(yIonIndexesTheo))},
                {Y_MZ_THEO, QVariant(qVectorToQByteArray(yIonMzValsTheo))},
                {Y_INTZ_THEO, QVariant(qVectorToQByteArray(yIonIntesitiesTheo))},
                {Y_IND_ACT, QVariant(qVectorToQByteArray(yIonIndexesActual))},
                {Y_MZ_ACT, QVariant(qVectorToQByteArray(yIonMzValsActual))},
                {Y_INTZ_ACT, QVariant(qVectorToQByteArray(yIonIntesitiesActual))},
                {Y_HILL_COS_SIM, QVariant(qVectorToQByteArray(yIonHillCosineSims))},

                {Y2_IND_THEO, QVariant(qVectorToQByteArray(y2IonIndexesTheo))},
                {Y2_MZ_THEO, QVariant(qVectorToQByteArray(y2IonMzValsTheo))},
                {Y2_INTZ_THEO, QVariant(qVectorToQByteArray(y2IonIntesitiesTheo))},
                {Y2_IND_ACT, QVariant(qVectorToQByteArray(y2IonIndexesActual))},
                {Y2_MZ_ACT, QVariant(qVectorToQByteArray(y2IonMzValsActual))},
                {Y2_INTZ_ACT, QVariant(qVectorToQByteArray(y2IonIntesitiesActual))},
                {Y2_HILL_COS_SIM, QVariant(qVectorToQByteArray(y2IonHillCosineSims))},

                {YNH3_IND_THEO, QVariant(qVectorToQByteArray(yNH3IonIndexesTheo))},
                {YNH3_MZ_THEO, QVariant(qVectorToQByteArray(yNH3IonMzValsTheo))},
                {YNH3_INTZ_THEO, QVariant(qVectorToQByteArray(yNH3IonIntesitiesTheo))},
                {YNH3_IND_ACT, QVariant(qVectorToQByteArray(yNH3IonIndexesActual))},
                {YNH3_MZ_ACT, QVariant(qVectorToQByteArray(yNH3IonMzValsActual))},
                {YNH3_INTZ_ACT, QVariant(qVectorToQByteArray(yNH3IonIntesitiesActual))},
                {YNH3_HILL_COS_SIM, QVariant(qVectorToQByteArray(yNH3IonHillCosineSims))},

                {YH2O_IND_THEO, QVariant(qVectorToQByteArray(yH2OIonIndexesTheo))},
                {YH2O_MZ_THEO, QVariant(qVectorToQByteArray(yH2OIonMzValsTheo))},
                {YH2O_INTZ_THEO, QVariant(qVectorToQByteArray(yH2OIonIntesitiesTheo))},
                {YH2O_IND_ACT, QVariant(qVectorToQByteArray(yH2OIonIndexesActual))},
                {YH2O_MZ_ACT, QVariant(qVectorToQByteArray(yH2OIonMzValsActual))},
                {YH2O_INTZ_ACT, QVariant(qVectorToQByteArray(yH2OIonIntesitiesActual))},
                {YH2O_HILL_COS_SIM, QVariant(qVectorToQByteArray(yH2OIonHillCosineSims))},

                {B_IND_THEO, QVariant(qVectorToQByteArray(bIonIndexesTheo))},
                {B_MZ_THEO, QVariant(qVectorToQByteArray(bIonMzValsTheo))},
                {B_INTZ_THEO, QVariant(qVectorToQByteArray(bIonIntesitiesTheo))},
                {B_IND_ACT, QVariant(qVectorToQByteArray(bIonIndexesActual))},
                {B_MZ_ACT, QVariant(qVectorToQByteArray(bIonMzValsActual))},
                {B_INTZ_ACT, QVariant(qVectorToQByteArray(bIonIntesitiesActual))},
                {B_HILL_COS_SIM, QVariant(qVectorToQByteArray(bIonHillCosineSims))},

                {B2_IND_THEO, QVariant(qVectorToQByteArray(b2IonIndexesTheo))},
                {B2_MZ_THEO, QVariant(qVectorToQByteArray(b2IonMzValsTheo))},
                {B2_INTZ_THEO, QVariant(qVectorToQByteArray(b2IonIntesitiesTheo))},
                {B2_IND_ACT, QVariant(qVectorToQByteArray(b2IonIndexesActual))},
                {B2_MZ_ACT, QVariant(qVectorToQByteArray(b2IonMzValsActual))},
                {B2_INTZ_ACT, QVariant(qVectorToQByteArray(b2IonIntesitiesActual))},
                {B2_HILL_COS_SIM, QVariant(qVectorToQByteArray(b2IonHillCosineSims))},

                {BNH3_IND_THEO, QVariant(qVectorToQByteArray(bNH3IonIndexesTheo))},
                {BNH3_MZ_THEO, QVariant(qVectorToQByteArray(bNH3IonMzValsTheo))},
                {BNH3_INTZ_THEO, QVariant(qVectorToQByteArray(bNH3IonIntesitiesTheo))},
                {BNH3_IND_ACT, QVariant(qVectorToQByteArray(bNH3IonIndexesActual))},
                {BNH3_MZ_ACT, QVariant(qVectorToQByteArray(bNH3IonMzValsActual))},
                {BNH3_INTZ_ACT, QVariant(qVectorToQByteArray(bNH3IonIntesitiesActual))},
                {BNH3_HILL_COS_SIM, QVariant(qVectorToQByteArray(bNH3IonHillCosineSims))},

                {BH2O_IND_THEO, QVariant(qVectorToQByteArray(bH2OIonIndexesTheo))},
                {BH2O_MZ_THEO, QVariant(qVectorToQByteArray(bH2OIonMzValsTheo))},
                {BH2O_INTZ_THEO, QVariant(qVectorToQByteArray(bH2OIonIntesitiesTheo))},
                {BH2O_IND_ACT, QVariant(qVectorToQByteArray(bH2OIonIndexesActual))},
                {BH2O_MZ_ACT, QVariant(qVectorToQByteArray(bH2OIonMzValsActual))},
                {BH2O_INTZ_ACT, QVariant(qVectorToQByteArray(bH2OIonIntesitiesActual))},
                {BH2O_HILL_COS_SIM, QVariant(qVectorToQByteArray(bH2OIonHillCosineSims))},

                {PREC_IND_THEO, QVariant(qVectorToQByteArray(precursorIonIndexesTheo))},
                {PREC_MZ_THEO, QVariant(qVectorToQByteArray(precursorIonMzValsTheo))},
                {PREC_INTZ_THEO, QVariant(qVectorToQByteArray(precursorIonIntesitiesTheo))},
                {PREC_IND_ACT, QVariant(qVectorToQByteArray(precursorIonIndexesActual))},
                {PREC_MZ_ACT, QVariant(qVectorToQByteArray(precursorIonMzValsActual))},
                {PREC_INTZ_ACT, QVariant(qVectorToQByteArray(precursorIonIntesitiesActual))},
                {PREC_HILL_COS_SIM,QVariant(qVectorToQByteArray(precursorIonIntesitiesActual))},

                {A_MZ_STD_ACT, QVariant(qVectorToQByteArray(aIonMzStdActual))},
                {A_HILL_LEN_ACT, QVariant(qVectorToQByteArray(aIonHillLengthActual))},
                {Y_MZ_STD_ACT, QVariant(qVectorToQByteArray(yIonMzStdActual))},
                {Y_HILL_LEN_ACT, QVariant(qVectorToQByteArray(yIonHillLengthActual))},
                {Y2_MZ_STD_ACT, QVariant(qVectorToQByteArray(y2IonMzStdActual))},
                {Y2_HILL_LEN_ACT, QVariant(qVectorToQByteArray(y2IonHillLengthActual))},
                {YNH3_MZ_STD_ACT, QVariant(qVectorToQByteArray(yNH3IonMzStdActual))},
                {YNH3_HILL_LEN_ACT, QVariant(qVectorToQByteArray(yNH3IonHillLengthActual))},
                {YH2O_MZ_STD_ACT, QVariant(qVectorToQByteArray(yH2OIonMzStdActual))},
                {YH2O_HILL_LEN_ACT, QVariant(qVectorToQByteArray(yH2OIonHillLengthActual))},
                {B_MZ_STD_ACT, QVariant(qVectorToQByteArray(bIonMzStdActual))},
                {B_HILL_LEN_ACT, QVariant(qVectorToQByteArray(bIonHillLengthActual))},
                {B2_MZ_STD_ACT, QVariant(qVectorToQByteArray(b2IonMzStdActual))},
                {B2_HILL_LEN_ACT, QVariant(qVectorToQByteArray(b2IonHillLengthActual))},
                {BNH3_MZ_STD_ACT, QVariant(qVectorToQByteArray(bNH3IonMzStdActual))},
                {BNH3_HILL_LEN_ACT, QVariant(qVectorToQByteArray(bNH3IonHillLengthActual))},
                {BH2O_MZ_STD_ACT, QVariant(qVectorToQByteArray(bH2OIonMzStdActual))},
                {BH2O_HILL_LEN_ACT, QVariant(qVectorToQByteArray(bH2OIonHillLengthActual))},
                {PREC_MZ_STD_ACT, QVariant(qVectorToQByteArray(precursorIonMzStdActual))},
                {PREC_HILL_LEN_ACT, QVariant(qVectorToQByteArray(precursorIonHillLengthActual))},

                {A_ISO_COS_ACT, QVariant(qVectorToQByteArray(aIonIsotopologueCosineSimActual))},
                {A_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(aIonIsotopologueChargeActual))},
                {A_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(aIonIsotopologueIntensityActual))},
                {Y_ISO_COS_ACT, QVariant(qVectorToQByteArray(yIonIsotopologueCosineSimActual))},
                {Y_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(yIonIsotopologueChargeActual))},
                {Y_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(yIonIsotopologueIntensityActual))},
                {Y2_ISO_COS_ACT, QVariant(qVectorToQByteArray(y2IonIsotopologueCosineSimActual))},
                {Y2_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(y2IonIsotopologueChargeActual))},
                {Y2_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(y2IonIsotopologueIntensityActual))},
                {YNH3_ISO_COS_ACT, QVariant(qVectorToQByteArray(yNH3IonIsotopologueCosineSimActual))},
                {YNH3_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(yNH3IonIsotopologueChargeActual))},
                {YNH3_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(yNH3IonIsotopologueIntensityActual))},
                {YH2O_ISO_COS_ACT, QVariant(qVectorToQByteArray(yH2OIonIsotopologueCosineSimActual))},
                {YH2O_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(yH2OIonIsotopologueChargeActual))},
                {YH2O_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(yH2OIonIsotopologueIntensityActual))},
                {B_ISO_COS_ACT, QVariant(qVectorToQByteArray(bIonIsotopologueCosineSimActual))},
                {B_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(bIonIsotopologueChargeActual))},
                {B_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(bIonIsotopologueIntensityActual))},
                {B2_ISO_COS_ACT, QVariant(qVectorToQByteArray(b2IonIsotopologueCosineSimActual))},
                {B2_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(b2IonIsotopologueChargeActual))},
                {B2_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(b2IonIsotopologueIntensityActual))},
                {BNH3_ISO_COS_ACT, QVariant(qVectorToQByteArray(bNH3IonIsotopologueCosineSimActual))},
                {BNH3_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(bNH3IonIsotopologueChargeActual))},
                {BNH3_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(bNH3IonIsotopologueIntensityActual))},
                {BH2O_ISO_COS_ACT, QVariant(qVectorToQByteArray(bH2OIonIsotopologueCosineSimActual))},
                {BH2O_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(bH2OIonIsotopologueChargeActual))},
                {BH2O_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(bH2OIonIsotopologueIntensityActual))},
                {PREC_ISO_COS_ACT, QVariant(qVectorToQByteArray(precursorIonIsotopologueCosineSimActual))},
                {PREC_ISO_CHRG_ACT, QVariant(qVectorToQByteArray(precursorIonIsotopologueChargeActual))},
                {PREC_ISO_INTZ_ACT, QVariant(qVectorToQByteArray(precursorIonIsotopologueIntensityActual))},
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace FrameExtractsReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideStringWithMods = dataMap.value(PEP_WITH_MODS).toString();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        scanNumberApex = dataMap.value(SCAN_NUMBER).toInt();
        frameIndexApex = dataMap.value(FRAME_INDEX).toInt();
        cosineSimSumWeighted = dataMap.value(COSINE_SUM_WEIGHTED).toDouble();
        aIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(A_IND_THEO).toByteArray());
        aIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(A_MZ_THEO).toByteArray());
        aIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(A_INTZ_THEO).toByteArray());
        aIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(A_IND_ACT).toByteArray());
        aIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(A_MZ_ACT).toByteArray());
        aIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(A_INTZ_ACT).toByteArray());
        aIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(A_HILL_COS_SIM).toByteArray());
        yIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(Y_IND_THEO).toByteArray());
        yIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(Y_MZ_THEO).toByteArray());
        yIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(Y_INTZ_THEO).toByteArray());
        yIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(Y_IND_ACT).toByteArray());
        yIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(Y_MZ_ACT).toByteArray());
        yIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(Y_INTZ_ACT).toByteArray());
        yIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(Y_HILL_COS_SIM).toByteArray());
        y2IonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(Y2_IND_THEO).toByteArray());
        y2IonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(Y2_MZ_THEO).toByteArray());
        y2IonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(Y2_INTZ_THEO).toByteArray());
        y2IonIndexesActual = bytesArrayToQVector<int>(dataMap.value(Y2_IND_ACT).toByteArray());
        y2IonMzValsActual = bytesArrayToQVector<double>(dataMap.value(Y2_MZ_ACT).toByteArray());
        y2IonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(Y2_INTZ_ACT).toByteArray());
        y2IonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(Y2_HILL_COS_SIM).toByteArray());
        yNH3IonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(YNH3_IND_THEO).toByteArray());
        yNH3IonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(YNH3_MZ_THEO).toByteArray());
        yNH3IonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(YNH3_INTZ_THEO).toByteArray());
        yNH3IonIndexesActual = bytesArrayToQVector<int>(dataMap.value(YNH3_IND_ACT).toByteArray());
        yNH3IonMzValsActual = bytesArrayToQVector<double>(dataMap.value(YNH3_MZ_ACT).toByteArray());
        yNH3IonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(YNH3_INTZ_ACT).toByteArray());
        yNH3IonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(YNH3_HILL_COS_SIM).toByteArray());
        yH2OIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(YH2O_IND_THEO).toByteArray());
        yH2OIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(YH2O_MZ_THEO).toByteArray());
        yH2OIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(YH2O_INTZ_THEO).toByteArray());
        yH2OIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(YH2O_IND_ACT).toByteArray());
        yH2OIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(YH2O_MZ_ACT).toByteArray());
        yH2OIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(YH2O_INTZ_ACT).toByteArray());
        yH2OIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(YH2O_HILL_COS_SIM).toByteArray());
        bIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(B_IND_THEO).toByteArray());
        bIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(B_MZ_THEO).toByteArray());
        bIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(B_INTZ_THEO).toByteArray());
        bIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(B_IND_ACT).toByteArray());
        bIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(B_MZ_ACT).toByteArray());
        bIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(B_INTZ_ACT).toByteArray());
        bIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(B_HILL_COS_SIM).toByteArray());
        b2IonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(B2_IND_THEO).toByteArray());
        b2IonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(B2_MZ_THEO).toByteArray());
        b2IonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(B2_INTZ_THEO).toByteArray());
        b2IonIndexesActual = bytesArrayToQVector<int>(dataMap.value(B2_IND_ACT).toByteArray());
        b2IonMzValsActual = bytesArrayToQVector<double>(dataMap.value(B2_MZ_ACT).toByteArray());
        b2IonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(B2_INTZ_ACT).toByteArray());
        b2IonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(B2_HILL_COS_SIM).toByteArray());
        bNH3IonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(BNH3_IND_THEO).toByteArray());
        bNH3IonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(BNH3_MZ_THEO).toByteArray());
        bNH3IonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(BNH3_INTZ_THEO).toByteArray());
        bNH3IonIndexesActual = bytesArrayToQVector<int>(dataMap.value(BNH3_IND_ACT).toByteArray());
        bNH3IonMzValsActual = bytesArrayToQVector<double>(dataMap.value(BNH3_MZ_ACT).toByteArray());
        bNH3IonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(BNH3_INTZ_ACT).toByteArray());
        bNH3IonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(BNH3_HILL_COS_SIM).toByteArray());
        bH2OIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(BH2O_IND_THEO).toByteArray());
        bH2OIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(BH2O_MZ_THEO).toByteArray());
        bH2OIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(BH2O_INTZ_THEO).toByteArray());
        bH2OIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(BH2O_IND_ACT).toByteArray());
        bH2OIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(BH2O_MZ_ACT).toByteArray());
        bH2OIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(BH2O_INTZ_ACT).toByteArray());
        bH2OIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(BNH3_HILL_COS_SIM).toByteArray());
        precursorIonIndexesTheo = bytesArrayToQVector<int>(dataMap.value(PREC_IND_THEO).toByteArray());
        precursorIonMzValsTheo = bytesArrayToQVector<double>(dataMap.value(PREC_MZ_THEO).toByteArray());
        precursorIonIntesitiesTheo = bytesArrayToQVector<double>(dataMap.value(PREC_INTZ_THEO).toByteArray());
        precursorIonIndexesActual = bytesArrayToQVector<int>(dataMap.value(PREC_IND_ACT).toByteArray());
        precursorIonMzValsActual = bytesArrayToQVector<double>(dataMap.value(PREC_MZ_ACT).toByteArray());
        precursorIonIntesitiesActual = bytesArrayToQVector<double>(dataMap.value(PREC_INTZ_ACT).toByteArray());
        precursorIonHillCosineSims = bytesArrayToQVector<double>(dataMap.value(PREC_HILL_COS_SIM).toByteArray());

        aIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(A_MZ_STD_ACT).toByteArray());
        aIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(A_HILL_LEN_ACT).toByteArray());
        yIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(Y_MZ_STD_ACT).toByteArray());
        yIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(Y_HILL_LEN_ACT).toByteArray());
        y2IonMzStdActual = bytesArrayToQVector<double>(dataMap.value(Y2_MZ_STD_ACT).toByteArray());
        y2IonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(Y2_HILL_LEN_ACT).toByteArray());
        yNH3IonMzStdActual = bytesArrayToQVector<double>(dataMap.value(YNH3_MZ_STD_ACT).toByteArray());
        yNH3IonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(YNH3_HILL_LEN_ACT).toByteArray());
        yH2OIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(YH2O_MZ_STD_ACT).toByteArray());
        yH2OIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(YH2O_HILL_LEN_ACT).toByteArray());
        bIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(B_MZ_STD_ACT).toByteArray());
        bIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(B_HILL_LEN_ACT).toByteArray());
        b2IonMzStdActual = bytesArrayToQVector<double>(dataMap.value(B2_MZ_STD_ACT).toByteArray());
        b2IonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(B2_HILL_LEN_ACT).toByteArray());
        bNH3IonMzStdActual = bytesArrayToQVector<double>(dataMap.value(BNH3_MZ_STD_ACT).toByteArray());
        bNH3IonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(BNH3_HILL_LEN_ACT).toByteArray());
        bH2OIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(BH2O_MZ_STD_ACT).toByteArray());
        bH2OIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(BH2O_HILL_LEN_ACT).toByteArray());
        precursorIonMzStdActual = bytesArrayToQVector<double>(dataMap.value(PREC_MZ_STD_ACT).toByteArray());
        precursorIonHillLengthActual = bytesArrayToQVector<int>(dataMap.value(PREC_HILL_LEN_ACT).toByteArray());

        aIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(A_ISO_COS_ACT).toByteArray());
        aIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(A_ISO_CHRG_ACT).toByteArray());
        aIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(A_ISO_INTZ_ACT).toByteArray());
        yIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(Y_ISO_COS_ACT).toByteArray());
        yIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(Y_ISO_CHRG_ACT).toByteArray());
        yIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(Y_ISO_INTZ_ACT).toByteArray());
        y2IonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(Y2_ISO_COS_ACT).toByteArray());
        y2IonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(Y2_ISO_CHRG_ACT).toByteArray());
        y2IonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(Y2_ISO_INTZ_ACT).toByteArray());
        yNH3IonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(YNH3_ISO_COS_ACT).toByteArray());
        yNH3IonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(YNH3_ISO_CHRG_ACT).toByteArray());
        yNH3IonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(YNH3_ISO_INTZ_ACT).toByteArray());
        yH2OIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(YH2O_ISO_COS_ACT).toByteArray());
        yH2OIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(YH2O_ISO_CHRG_ACT).toByteArray());
        yH2OIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(YH2O_ISO_INTZ_ACT).toByteArray());
        bIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(B_ISO_COS_ACT).toByteArray());
        bIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(B_ISO_CHRG_ACT).toByteArray());
        bIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(B_ISO_INTZ_ACT).toByteArray());
        b2IonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(B2_ISO_COS_ACT).toByteArray());
        b2IonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(B2_ISO_CHRG_ACT).toByteArray());
        b2IonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(B2_ISO_INTZ_ACT).toByteArray());
        bNH3IonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(BNH3_ISO_COS_ACT).toByteArray());
        bNH3IonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(BNH3_ISO_CHRG_ACT).toByteArray());
        bNH3IonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(BNH3_ISO_INTZ_ACT).toByteArray());
        bH2OIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(BH2O_ISO_COS_ACT).toByteArray());
        bH2OIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(BH2O_ISO_CHRG_ACT).toByteArray());
        bH2OIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(BH2O_ISO_INTZ_ACT).toByteArray());
        precursorIonIsotopologueCosineSimActual = bytesArrayToQVector<double>(dataMap.value(PREC_ISO_COS_ACT).toByteArray());
        precursorIonIsotopologueChargeActual = bytesArrayToQVector<int>(dataMap.value(PREC_ISO_CHRG_ACT).toByteArray());
        precursorIonIsotopologueIntensityActual = bytesArrayToQVector<double>(dataMap.value(PREC_ISO_INTZ_ACT).toByteArray());

        ERR_RETURN
    }

};


#endif //PYTHIADIACPP_FRAMEEXTRACTSREADER_H
