//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePair.h"

#include "BiophysicalCalcs.h"
#include "MsUtils.h"

TargetDecoyCandidatePair::TargetDecoyCandidatePair()
: m_peptideString("")
, m_decoyMassDelta(-1.0)
, m_fragLibReaderRowPntr(nullptr)
, m_decoySharesSequenceWithOtherTarget(false)
{}

Err TargetDecoyCandidatePair::init() {

	ERR_INIT

	e = ErrorUtils::isFalse(m_fragLibReaderRowPntr == nullptr); ree;

	int chargeNotUsed;
	e = MsUtils::peptideStringWithModsFromPeptideSequenceChargeKey(
		m_fragLibReaderRowPntr->peptideSequenceChargeKey,
		&m_peptideStringWithMods,
		&chargeNotUsed
		); ree;

	m_peptideString = m_peptideStringWithMods.removeUniModChars();

	e = buildResidueMutations(); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePair::buildResidueMutations() {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_peptideString); ree;

	const QMap<QChar, double> diannMutateAminoAcidToMass = AminoAcids::diannMutateAminoAcidToMass();
	const QMap<QChar, QChar> diannMutateAminoAcidToAminoAcid = AminoAcids::diannMutateAminoAcidToResidue();

	constexpr int firstMutationResidue = 1;
	ResidueMutation resMutFirst;
	resMutFirst.index = firstMutationResidue;
	resMutFirst.residueMutatedTo = diannMutateAminoAcidToAminoAcid.value(m_peptideString.at(firstMutationResidue));
	resMutFirst.residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(firstMutationResidue));

	const int penultimateResidue =  m_peptideString.size() - 2;
	ResidueMutation resMutPen;
	resMutPen.index = penultimateResidue;
	resMutPen.residueMutatedTo = diannMutateAminoAcidToAminoAcid.value(m_peptideString.at(penultimateResidue));
	resMutPen.residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(penultimateResidue));

	// m_residueMutations.clear();
	// m_residueMutations = {
	// 	resMutFirst,
	// 	resMutPen
	// };

	constexpr float minAbsoluteDeltaMass = 15.0;

	m_decoyMassDelta = calculateDecoyMass();

	// if (m_decoyMassDelta < minAbsoluteDeltaMass) {
		constexpr int firstMutationResidue2 = 2;
		ResidueMutation resMutFirst2;
		resMutFirst2.index = firstMutationResidue2;
		resMutFirst2.residueMutatedTo = diannMutateAminoAcidToAminoAcid.value(m_peptideString.at(firstMutationResidue2));
		resMutFirst2.residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(firstMutationResidue2));

		const int penultimateResidue2 =  m_peptideString.size() - 3;
		ResidueMutation resMutPen2;
		resMutPen2.index = penultimateResidue2;
		resMutPen2.residueMutatedTo = diannMutateAminoAcidToAminoAcid.value(m_peptideString.at(penultimateResidue2));
		resMutPen2.residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(penultimateResidue2));
		m_residueMutations.clear();
		m_residueMutations.append({
			resMutFirst,
			// resMutPen,
			// resMutFirst2,
			resMutPen2
		});
		m_decoyMassDelta = calculateDecoyMass();
	//
	// 	float deltaMassBest  = m_decoyMassDelta;
	// 	int bestResidueIndex = -1;
	// 	for (int i = firstMutationResidue; i < penultimateResidue; ++i) {
	// 		const double residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(i));
	// 		const float deltaMass = std::abs(m_decoyMassDelta + residueDeltaMass);
	//
	// 		if (deltaMass > deltaMassBest) {
	// 			deltaMassBest = deltaMass;
	// 			bestResidueIndex = i;
	// 		}
	//
	// 		if (deltaMassBest > minAbsoluteDeltaMass) {
	// 			break;
	// 		}
	// 	}
	//
	// 	ResidueMutation resMut;
	// 	resMut.index = bestResidueIndex;
	// 	resMut.residueMutatedTo = diannMutateAminoAcidToAminoAcid.value(m_peptideString.at(bestResidueIndex));
	// 	resMut.residueDeltaMass = diannMutateAminoAcidToMass.value(m_peptideString.at(bestResidueIndex));
	// 	m_residueMutations.push_back(resMut);
	//
	// 	m_decoyMassDelta = calculateDecoyMass();
	// }

	ERR_RETURN
}

float TargetDecoyCandidatePair::calculateDecoyMass() const {
	return std::accumulate(
		m_residueMutations.begin(),
		m_residueMutations.end(),
		0.0f,
		[](float sum, const ResidueMutation &rm){return sum + rm.residueDeltaMass;});
}

void TargetDecoyCandidatePair::setFragLibReaderRowPntr(FragLibReaderRow *fragLibReaderRowPntr) {
    m_fragLibReaderRowPntr = fragLibReaderRowPntr;
}

QString TargetDecoyCandidatePair::proteinGroups() const {
    return m_fragLibReaderRowPntr->proteinGroups;
}

PeptideString TargetDecoyCandidatePair::peptideString() const {
    return m_peptideString;
}

PeptideString TargetDecoyCandidatePair::peptideStringDecoy() const {
	PeptideString peptideStringDecoy = m_peptideString;
	for (const ResidueMutation &rm : m_residueMutations) {
		peptideStringDecoy[rm.index] = rm.residueMutatedTo;
	}

	return peptideStringDecoy;
}

PeptideStringWithMods TargetDecoyCandidatePair::peptideStringWithMods() const {
    return m_peptideStringWithMods;
}

PeptideStringWithMods TargetDecoyCandidatePair::peptideStringWithModsDecoy() const {
	PeptideStringWithMods peptideStringWithMods = m_peptideStringWithMods;


	for (const ResidueMutation &rm : m_residueMutations) {

		bool isUiModOn = false;
		int currentResidueIndex = 0;

		for (int currentResidueWithModsIndex = 0; currentResidueWithModsIndex < peptideStringWithMods.size(); ++currentResidueWithModsIndex) {

			const QChar &ch = peptideStringWithMods[currentResidueWithModsIndex];

			if (ch == '[' || ch == '(') {
				isUiModOn = true;
			}
			if (ch == ']' || ch == ')') {
				isUiModOn = false;
				continue;
			}

			if (isUiModOn) {
				continue;
			}

			if (currentResidueIndex == rm.index) {
				peptideStringWithMods[currentResidueWithModsIndex] = rm.residueMutatedTo;
				break;
			}

			currentResidueIndex++;

		}
	}

	return peptideStringWithMods;
}

namespace {

    QVector<MS2Ion> buildMS2Ions(const FragLibReaderRow *flrr) {

        const QStringList ionLabelsSplit = flrr->ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, Qt::SkipEmptyParts);

        QVector<MS2Ion> ms2IonsBuilder;
        ms2IonsBuilder.reserve(flrr->mzVals.size());
        for (int i = 0; i < flrr->mzVals.size(); i++) {
            MS2Ion ms2Ion;
            ms2Ion.mz = flrr->mzVals.at(i);
            ms2Ion.intensity = flrr->intensityVals.at(i);
            ms2Ion.ionLabel = ionLabelsSplit.at(i);
            ms2Ion.charge = ms2Ion.ionLabel.contains("^2") ? 2 : 1;

            ms2IonsBuilder.push_back(ms2Ion);
        }

        constexpr float mzMin = 200.0;
        constexpr float mzMax = 1500.0;
        MS2Ion::filterMS2IonsByMz(
                mzMin,
                mzMax,
                &ms2IonsBuilder
                );

        MS2Ion::sortMS2IonsIntensityDesc(&ms2IonsBuilder);
        for (int i = 0; i < ms2IonsBuilder.size(); i++) {
            MS2Ion &ms2Ion = ms2IonsBuilder[i];
            ms2Ion.rank = i;
        }

        if (ms2IonsBuilder.isEmpty()) {
            qDebug()
            << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "No MS2 ions found for"
            << flrr->peptideSequenceChargeKey
            << flrr->mzVals
            << flrr->intensityVals
            << flrr->ionLabels;
        }

        return ms2IonsBuilder;
    }

}//namespace
QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsTarget() {

	if (m_ms2Targets.isEmpty()) {
		m_ms2Targets = buildMS2Ions(m_fragLibReaderRowPntr);
	}

    return m_ms2Targets;
}

QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsDecoy() {

    QVector<MS2Ion> ms2IonsDec = mutateCandidatePeptideTarget();

    return ms2IonsDec;
}

QVector<MS2Ion> TargetDecoyCandidatePair::mutateCandidatePeptideTarget() {

    ERR_INIT

	e = ErrorUtils::isNotEmpty(m_residueMutations); rree;

    QVector<MS2Ion> ms2IonDecoys = ms2IonsTarget();
	const int pepLength = m_peptideString.size();

	for (const ResidueMutation &rm : m_residueMutations) {

		for (MS2Ion &ms2Ion : ms2IonDecoys) {

			QPair<IonIndex, IonType> ionLableInfo;
			e = ms2Ion.getIonLabelInfo(&ionLableInfo);

			if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

				if (ionLableInfo.second.contains("^2")) {

					if (ionLableInfo.first - 1  >= ionLableInfo.second.size()) {
						ms2Ion.mz += rm.residueDeltaMass / 2.0f;

					}
				}
				else {

					if (ionLableInfo.first - 1  >= rm.index) {
						ms2Ion.mz += rm.residueDeltaMass;
					}
				}
			}

			else if (ionLableInfo.second.contains('y')) {

				const int bToYIndex = pepLength - rm.index - 1;

				if (ionLableInfo.second.contains("^2")) {

					if (ionLableInfo.first - 1 >= bToYIndex) {
						ms2Ion.mz += rm.residueDeltaMass / 2.0f; //NOTE: do not static cast to float
					}

				}
				else {
					if (ionLableInfo.first - 1 >= bToYIndex) {
						ms2Ion.mz += rm.residueDeltaMass; //NOTE: do not static cast to float
					}

				}
			}

			else {
				qDebug() << "Non b/y/a ion" << ionLableInfo;
			}

		}
	}

    return ms2IonDecoys;
}

float TargetDecoyCandidatePair::mz(bool isDecoy) const {
    return BiophysicalCalcs::calculateThomsonFromMass(mass(isDecoy) , charge());
}

bool TargetDecoyCandidatePair::isDecoy() const {
	return m_fragLibReaderRowPntr->isDecoy;
}

int TargetDecoyCandidatePair::charge() const {
    return m_fragLibReaderRowPntr->precursorCharge;
}

float TargetDecoyCandidatePair::mass(bool isDecoy) const {
    return isDecoy ? static_cast<float>(m_fragLibReaderRowPntr->mass + m_decoyMassDelta)
					: static_cast<float>(m_fragLibReaderRowPntr->mass);
}

float TargetDecoyCandidatePair::iRt(bool isDecoy) const {

	float decoyAdjustment = 0;
	if (isDecoy) {
		const PeptideString ps = peptideString();
		const QChar secondAA = ps[1];
		const QChar penultimateAA = ps[ps.size() - 2];

		decoyAdjustment += UniModNamespace::iRtAdjustments.value(secondAA)
						+ UniModNamespace::iRtAdjustments.value(penultimateAA);
	}

    return static_cast<float>(m_fragLibReaderRowPntr->iRT) + decoyAdjustment;

}

float TargetDecoyCandidatePair::iIM() const {
    return static_cast<float>(m_fragLibReaderRowPntr->iM);
}

bool TargetDecoyCandidatePair::isInit() const {
	return m_fragLibReaderRowPntr != nullptr && !m_peptideString.isEmpty();
}

int TargetDecoyCandidatePair::totalFragmentCount() const {
    return m_fragLibReaderRowPntr->mzVals.size();
}

void TargetDecoyCandidatePair::mutateCandidatePeptideTargetTestAccess(
	const PeptideStringWithMods &peptideStringWithMods,
	const QVector<MS2Ion> &ms2IonTarget
	) {
	//
	// qDebug() << peptideStringWithMods.bSeries(1, AminoAcids());
	// qDebug() << peptideStringWithMods.ySeries(1, AminoAcids());
	//
	// const QVector<MS2Ion> decoys = mutateCandidatePeptideTarget();
	// qDebug() << ms2IonTarget;
	// qDebug() << decoys;
}

void TargetDecoyCandidatePair::mangleMs2IonsDecoy(QVector<MS2Ion> *ms2Ions) {
    for (MS2Ion &ms2Ion : *ms2Ions) {
        ms2Ion.mz += 0.1;
    }
}

void TargetDecoyCandidatePair::decoySharesSequenceWithOtherTarget(bool val) {
    m_decoySharesSequenceWithOtherTarget = val;
}
