//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePair.h"

#include "BiophysicalCalcs.h"

TargetDecoyCandidatePair::TargetDecoyCandidatePair()
: m_peptideStringWithMods("")
, m_decoyMassDelta(-1.0)
, m_fragLibReaderRowPntr(nullptr)
, m_decoySharesSequenceWithOtherTarget(false)
{}

TargetDecoyCandidatePair::TargetDecoyCandidatePair(
        const PeptideStringWithMods &peptideStringWithMods,
        float decoyMassDelta
)
: m_peptideStringWithMods(peptideStringWithMods)
, m_decoyMassDelta(decoyMassDelta)
, m_fragLibReaderRowPntr(nullptr)
, m_decoySharesSequenceWithOtherTarget(false)
{}

void TargetDecoyCandidatePair::setFragLibReaderRowPntr(FragLibReaderRow *fragLibReaderRowPntr) {
    m_fragLibReaderRowPntr = fragLibReaderRowPntr;
}

QString TargetDecoyCandidatePair::proteinGroups() const {
    return m_fragLibReaderRowPntr->proteinGroups;
}

PeptideString TargetDecoyCandidatePair::peptideString() const {
    return m_peptideStringWithMods.removeUniModChars();
}

PeptideStringWithMods TargetDecoyCandidatePair::peptideStringWithMods() const {
    return m_peptideStringWithMods;
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
QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsTarget() const {
    return buildMS2Ions(m_fragLibReaderRowPntr);
}

namespace {

    QVector<MS2Ion> mutateCandidatePeptideTarget(
        const PeptideStringWithMods &peptideStringWithMods,
        const QVector<MS2Ion> &ms2IonTarget
        ) {

        ERR_INIT

        const QMap<QChar, double> diannMutateAminoAcidToMass = AminoAcids::diannMutateAminoAcidToMass();

        const PeptideString &peptideString = peptideStringWithMods.removeUniModChars();

        constexpr int firstIndexToMutate = 1;
        const int secondIndexToMutate = peptideString.size() - 2;

        const double nTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(firstIndexToMutate));
        const double cTermDeltaMass = diannMutateAminoAcidToMass.value(peptideString.at(secondIndexToMutate));
        const double nTermDeltaMassCharge2 = nTermDeltaMass / 2.0;
        const double cTermDeltaMassCharge2 = cTermDeltaMass / 2.0;

        QVector<MS2Ion> ms2IonDecoys;

        for (const MS2Ion &ms2Ion : ms2IonTarget) {

            MS2Ion ms2IonDecoy = ms2Ion;

            QPair<IonIndex, IonType> ionLableInfo;
            e = ms2IonDecoy.getIonLabelInfo(&ionLableInfo);

            if (ionLableInfo.second.contains('b') || ionLableInfo.second.contains('a')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first - 1  >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2;
                    }

                    if (ionLableInfo.first - 1  >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2;
                    }
                }
                else {

                    if (ionLableInfo.first - 1  >= firstIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass;
                    }

                    if (ionLableInfo.first - 1  >= secondIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass;
                    }
                }
            }

            else if (ionLableInfo.second.contains('y')) {

                if (ionLableInfo.second.contains("^2")) {

                    if (ionLableInfo.first - 1  >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first - 1  >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMassCharge2; //NOTE: do not static cast to float
                    }
                }
                else {

                    if (ionLableInfo.first - 1  >= firstIndexToMutate) {
                        ms2IonDecoy.mz += cTermDeltaMass; //NOTE: do not static cast to float
                    }

                    if (ionLableInfo.first - 1 >= secondIndexToMutate) {
                        ms2IonDecoy.mz += nTermDeltaMass; //NOTE: do not static cast to float
                    }
                }
            }

            else {
                qDebug() << "Non b/y/a ion" << ionLableInfo;
            }

            ms2IonDecoys.push_back(ms2IonDecoy);
        }

        return ms2IonDecoys;
    }

}//namespace
QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsDecoy() const {

    QVector<MS2Ion> ms2IonsDec = mutateCandidatePeptideTarget(peptideStringWithMods(), ms2IonsTarget());
    if (m_decoySharesSequenceWithOtherTarget) {
        mangleMs2IonsDecoy(&ms2IonsDec);
    }
    return ms2IonsDec;
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

int TargetDecoyCandidatePair::totalFragmentCount() const {
    return m_fragLibReaderRowPntr->mzVals.size();
}

void TargetDecoyCandidatePair::setPeptideStringWithMods(const PeptideStringWithMods &peptideStringWithMods) {
    m_peptideStringWithMods = peptideStringWithMods;
}

void TargetDecoyCandidatePair::mutateCandidatePeptideTargetTestAccess(
	const PeptideStringWithMods &peptideStringWithMods,
	const QVector<MS2Ion> &ms2IonTarget
	) {

	qDebug() << peptideStringWithMods.bSeries(1, AminoAcids());
	qDebug() << peptideStringWithMods.ySeries(1, AminoAcids());

	const QVector<MS2Ion> decoys = mutateCandidatePeptideTarget(peptideStringWithMods, ms2IonTarget);
	qDebug() << ms2IonTarget;
	qDebug() << decoys;
}

void TargetDecoyCandidatePair::mangleMs2IonsDecoy(QVector<MS2Ion> *ms2Ions) {
    for (MS2Ion &ms2Ion : *ms2Ions) {
        ms2Ion.mz += 0.1;
    }
}

void TargetDecoyCandidatePair::decoySharesSequenceWithOtherTarget(bool val) {
    m_decoySharesSequenceWithOtherTarget = val;
}


