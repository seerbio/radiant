//
// Created by andrewnichols on 9/5/25.
//

#include "PrecursorMzArrangetron.h"

#include "ErrorUtils.h"

Err PrecursorMzArrangetron::init(MsReaderPointerAcc *msReaderPtr) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;
	m_msReaderPtr = msReaderPtr;

	e = extractScansParallel(); ree;


	ERR_RETURN
}

namespace {

	constexpr int precursorMzKeyHashingPrecision = 3;

	void sortScanNumberMzIntensitiesMzPrecursorMzAscMS2Asc(QVector<ScanNumberMzIntensity> *scanNumberMzIntensities) {

		const auto sortLogicMzPrecursorMzAscMS2Asc = [](const ScanNumberMzIntensity &l, const ScanNumberMzIntensity &r) {
			if (MathUtils::hashDecimal(l.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision)
				== MathUtils::hashDecimal(r.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision)
				) {
				return l.mzVal < r.mzVal;
			}
			return l.scanInfoPntr->precursorTargetMz < r.scanInfoPntr->precursorTargetMz;
		};

		std::sort(
			scanNumberMzIntensities->begin(),
			scanNumberMzIntensities->end(),
			sortLogicMzPrecursorMzAscMS2Asc
			);
	}

	void sortMzMS2Asc(QVector<ScanNumberMzIntensity> &v) {
		std::sort(
			v.begin(),
			v.end(),
			[](const ScanNumberMzIntensity &l, const ScanNumberMzIntensity &r) {
				return l.mzVal < r.mzVal;
			});
	}

	QPair<Err, QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>>> extractScansParallelLogic(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		) {
		ERR_INIT

		QVector<MsScan> msScans;
		e = msReaderPtr->ptr->extractScanPoints(
			scanInfosPntrs,
			&msScans
			); rree;

		const int tranchePointCount = std::accumulate(
			msScans.begin(),
			msScans.end(),
			0,
			[](int sum, const MsScan &msScan){return sum + msScan.msScanInfoPntr->pointCount;}
			);

		QVector<ScanNumberMzIntensity> scanNumberMzIntensities;
		scanNumberMzIntensities.reserve(tranchePointCount);

		for (const MsScan &scan : msScans) {

			for (int i = 0; i < scan.msScanInfoPntr->pointCount; i++) {

				constexpr float mzMS2Min = 300; //TODO get values from MsReader for these vals
				constexpr float mzMS2Max = 1600; //TODO get values from MsReader for these vals
				const float mzMS2Val = scan.mzVals[i];

				if (mzMS2Val < mzMS2Min || mzMS2Val > mzMS2Max) {
					continue;
				}

				ScanNumberMzIntensity sni;
				sni.scanInfoPntr = scan.msScanInfoPntr;
				sni.mzVal = mzMS2Val;
				sni.intensityVal = scan.intensityVals[i];
				scanNumberMzIntensities.push_back(sni);
			}
		}
		sortScanNumberMzIntensitiesMzPrecursorMzAscMS2Asc(&scanNumberMzIntensities);

		QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>> precursorMzKeyVsScanNumberMzIntensities;
		QVector<ScanNumberMzIntensity> scanNumberMzIntensityLoader;
		constexpr int reserveEstimation = 1e3;
		scanNumberMzIntensityLoader.reserve(reserveEstimation);

		int currentMzPrecursorHashed = -1;
		for (const ScanNumberMzIntensity &snmi : scanNumberMzIntensities) {

			if (MathUtils::hashDecimal(snmi.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision)
				!= currentMzPrecursorHashed) {

				if (!scanNumberMzIntensityLoader.isEmpty()) {

					sortMzMS2Asc(scanNumberMzIntensityLoader);

					precursorMzKeyVsScanNumberMzIntensities[currentMzPrecursorHashed].append(scanNumberMzIntensityLoader);
					scanNumberMzIntensityLoader.clear();
				}

				currentMzPrecursorHashed = MathUtils::hashDecimal(snmi.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision);
			}

			scanNumberMzIntensityLoader.push_back(snmi);
		}

		if (!precursorMzKeyVsScanNumberMzIntensities.isEmpty()) {
			const PrecursorMzKey precursorKey = MathUtils::hashDecimal(
				scanNumberMzIntensityLoader.front().scanInfoPntr->precursorTargetMz,
				precursorMzKeyHashingPrecision
				);
			precursorMzKeyVsScanNumberMzIntensities[precursorKey].append(scanNumberMzIntensityLoader);
		}

		return {e, precursorMzKeyVsScanNumberMzIntensities};
	}

}//namespace
Err PrecursorMzArrangetron::extractScansParallel() {

	ERR_INIT

	e = ErrorUtils::isTrue(m_msReaderPtr->isInit()); ree;

	constexpr int msLevel = 2;
	const QMap<ScanNumber, MsScanInfo*> scanNumberVsMsScanInfoPntrs = m_msReaderPtr->ptr->getMsScanInfos(msLevel);

	QVector<MsScanInfo*> msScanInfosPntrsSortedPrecursorAsc = scanNumberVsMsScanInfoPntrs.values().toVector();
	std::sort(
		msScanInfosPntrsSortedPrecursorAsc.begin(),
		msScanInfosPntrsSortedPrecursorAsc.end(),
		[](const MsScanInfo *l, const MsScanInfo *r){return l->precursorTargetMz < r->precursorTargetMz;}
		);

	float lp = -1;
	for (const MsScanInfo* msi : msScanInfosPntrsSortedPrecursorAsc) {

		if (lp > msi->precursorTargetMz) {
			qDebug() << lp << msi->precursorTargetMz << "SDLFJS";
		}

		lp = msi->precursorTargetMz;
	}

	constexpr int threads = 128;
	QVector<QVector<MsScanInfo*>> scanInfosPntrsTranched;
	e = ParallelUtils::trancheVectorForParallelizationInOrder(
		msScanInfosPntrsSortedPrecursorAsc,
		threads,
		0,
		&scanInfosPntrsTranched
		); ree;

#define PARALLEL_LOAD_666
#ifdef PARALLEL_LOAD_666
	const auto binderLogic = std::bind(
		extractScansParallelLogic,
		std::placeholders::_1,
		m_msReaderPtr
		);

	QFuture<QPair<Err, QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>>>> future = QtConcurrent::mapped(
		scanInfosPntrsTranched,
		binderLogic
		);
	future.waitForFinished();

	for (const QPair<Err, QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>>> &res : future) {
		e = res.first; ree;
		const QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>> m = res.second;
		for (auto it = m.begin(); it != m.end(); it++) {
			m_precursorMzKeyVsScanNumberMzIntensity[it.key()].append(it.value());
		}
	}

	QFuture<void> future2 = QtConcurrent::map(
		m_precursorMzKeyVsScanNumberMzIntensity,
		sortMzMS2Asc
		);
	future2.waitForFinished();

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsPoints extracted";
	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Precursor target count" << m_precursorMzKeyVsScanNumberMzIntensity.size();
#else
	for (const QVector<MsScanInfo*> &msi : scanInfosPntrsTranched) {
		const QPair<Err, QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>>> res = extractScansParallelLogic(msi, m_msReaderPtr);
	}
#endif

// #define SANITY_CHECK
#ifdef SANITY_CHECK
	float lastMzPrecursor = -1;
	float lastMzMS2 = -1;
	for (auto it = m_precursorMzKeyVsScanNumberMzIntensity.begin(); it != m_precursorMzKeyVsScanNumberMzIntensity.end(); ++it) {

		int key = it.key();
		const QVector<ScanNumberMzIntensity> &v = it.value();

		// qDebug()
		// << qPrintable(S_GLOBAL_TIMER.elapsed())
		// << "Precursor Key" << key
		// << "Precursor Mz" << v.front().scanInfoPntr->precursorTargetMz
		// << "Point count" << v.size();

		for (const ScanNumberMzIntensity &snmi : v) {
			const float currentMzPrecursor = snmi.scanInfoPntr->precursorTargetMz;
			const float currentMzMS2 = snmi.mzVal;

			e = ErrorUtils::isTrue(
				MathUtils::hashDecimal(lastMzPrecursor, precursorMzKeyHashingPrecision) == MathUtils::hashDecimal(currentMzPrecursor, precursorMzKeyHashingPrecision)
				|| lastMzPrecursor < currentMzPrecursor
				);
			if (e != eNoError) {
				qDebug() << lastMzPrecursor << currentMzPrecursor;
				rrr(eValueError);
			}

			lastMzMS2 = MathUtils::hashDecimal(lastMzPrecursor, precursorMzKeyHashingPrecision)
						!= MathUtils::hashDecimal(currentMzPrecursor, precursorMzKeyHashingPrecision)
					  ? -1
					  : lastMzMS2;

			e = ErrorUtils::isTrue(
			MathUtils::hashDecimal(lastMzMS2, precursorMzKeyHashingPrecision)
					== MathUtils::hashDecimal(currentMzMS2, precursorMzKeyHashingPrecision) || lastMzMS2 < currentMzMS2
				);
			if (e != eNoError) {
				qDebug() << lastMzMS2 << currentMzMS2;
				rrr(eValueError);
			}

			lastMzPrecursor = currentMzPrecursor;
			lastMzMS2 = currentMzMS2;
		}
	}
#endif

	ERR_RETURN
}

