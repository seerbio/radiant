//
// Created by Drucifer on 3/21/2022.
//

#include "GlobalSettings.h"

const GlobalSettings S_GLOBAL_SETTINGS;

const GlobalTimer S_GLOBAL_TIMER;

QString GlobalSettings::VERSION() {
    return QStringLiteral("1.0");
}

namespace UniModNamespace {

	const QHash<QChar, float> iRtAdjustments = {
		{'G', 0.5},
		{'A', 0.5},
		{'V', 0.25},
		{'L', -0.25},
		{'X', 0.5},
		{'I', -0.25},
		{'F', 0.25},
		{'M', 0.25},
		{'P', 0.5},
		{'W', 0.25},
		{'S', 0.25},
		{'C', -0.25},
		{'U', -0.5},
		{'T', -0.25},
		{'Y', -0.5},
		{'H', -0.5},
		{'K', 0.5},
		{'R', 0.5},
		{'Q', -0.25},
		{'E', -0.25},
		{'N', 0.25},
		{'D', 0.25}
	};

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
}
