//
// Created by andrewnichols on 9/5/25.
//

#include "PrecursorMzArrangetron.h"

#include "ErrorUtils.h"

PrecursorMzArrangetron::PrecursorMzArrangetron() : m_msReaderPtr(nullptr) {}

Err PrecursorMzArrangetron::init(MsReaderPointerAcc *msReaderPtr) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;
	m_msReaderPtr = msReaderPtr;

	e = extractScansParallel(); ree;

	ERR_RETURN
}

bool PrecursorMzArrangetron::isInit() const {
	return !m_precursorMzKeyVsScanNumberMzIntensity.isEmpty();
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

	Err extractScansParallelLogic(
	    const QVector<MsScanInfo *> &scanInfosPntrs,
	    MsReaderPointerAcc *msReaderPtr,
	    QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>> *precursorMzKeyVsScanNumberMzIntensities,
	    QMutex *mutex // New parameter for thread safety
	) {
	    ERR_INIT

	    QVector<MsScan> msScans;
	    e = msReaderPtr->ptr->extractScanPoints(
	        scanInfosPntrs,
	        &msScans
	    );
	    ree;

	    const int tranchePointCount = std::accumulate(
	        msScans.begin(),
	        msScans.end(),
	        0,
	        [](int sum, const MsScan &msScan) { return sum + msScan.msScanInfoPntr->pointCount; }
	    );

	    QVector<ScanNumberMzIntensity> scanNumberMzIntensities;
	    scanNumberMzIntensities.reserve(tranchePointCount);

	    for (const MsScan &scan : msScans) {
	        for (int i = 0; i < scan.msScanInfoPntr->pointCount; i++) {
	            constexpr float mzMS2Min = 300;  // TODO: get values from MsReader
	            constexpr float mzMS2Max = 1600; // TODO: get values from MsReader
	            const float mzMS2Val = scan.mzVals[i];

	            // Skip invalid values
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

	    QVector<ScanNumberMzIntensity> scanNumberMzIntensityLoader;
	    constexpr int reserveEstimation = 1e3;
	    scanNumberMzIntensityLoader.reserve(reserveEstimation);

	    int currentMzPrecursorHashed = -1;
	    for (const ScanNumberMzIntensity &snmi : scanNumberMzIntensities) {
	        if (MathUtils::hashDecimal(snmi.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision)
	            != currentMzPrecursorHashed) {
	            if (!scanNumberMzIntensityLoader.isEmpty()) {
	                sortMzMS2Asc(scanNumberMzIntensityLoader);

	                // Lock the mutex while modifying the shared QMap
	                {
	                    QMutexLocker locker(mutex);
	                    (*precursorMzKeyVsScanNumberMzIntensities)[currentMzPrecursorHashed].append(scanNumberMzIntensityLoader);
	                }

	                scanNumberMzIntensityLoader.clear();
	            }

	            currentMzPrecursorHashed = MathUtils::hashDecimal(snmi.scanInfoPntr->precursorTargetMz, precursorMzKeyHashingPrecision);
	        }

	        scanNumberMzIntensityLoader.push_back(snmi);
	    }

	    if (!scanNumberMzIntensityLoader.isEmpty()) {
	        const PrecursorMzKey precursorKey = MathUtils::hashDecimal(
	            scanNumberMzIntensityLoader.front().scanInfoPntr->precursorTargetMz,
	            precursorMzKeyHashingPrecision
	        );

	        // Lock the mutex while modifying the shared QMap
	        {
	            QMutexLocker locker(mutex);
	            (*precursorMzKeyVsScanNumberMzIntensities)[precursorKey].append(scanNumberMzIntensityLoader);
	        }
	    }

	    ERR_RETURN
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
	QMutex mutex;

	const auto binderLogic = std::bind(
		extractScansParallelLogic,
		std::placeholders::_1,
		m_msReaderPtr,
		&m_precursorMzKeyVsScanNumberMzIntensity, // Shared data structure
		&mutex                                    // Mutex for thread safety
	);

	QFuture<Err> future = QtConcurrent::mapped(
		scanInfosPntrsTranched,
		binderLogic
	);
	future.waitForFinished();

	for (const Err res : future) {
		e = res; ree;
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

	ERR_RETURN
}

namespace {

	QVector<int> calculateTranchSize(int size, int threadCount) {

		const int tranchSize = size / threadCount;
		int trancheSizeRemainder = size % threadCount;
		QVector<int> trancheSizes(threadCount, tranchSize);

		int counter = 0;
		while (trancheSizeRemainder-- > 0) {
			trancheSizes[counter++]++;
		}

		return trancheSizes;
	}

}//namespace
Err PrecursorMzArrangetron::trancheMsScanPoints(
	int threadCount,
	QVector<QVector<ScanNumberMzIntensity*>> *scanNumberMzIntensityTranched
	) {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(threadCount); ree;
	e = ErrorUtils::isNotEmpty(m_precursorMzKeyVsScanNumberMzIntensity); ree;

	const int pointCountAllScans = std::accumulate(
		m_precursorMzKeyVsScanNumberMzIntensity.begin(),
		m_precursorMzKeyVsScanNumberMzIntensity.end(),
		0,
		[](int sum, const QVector<ScanNumberMzIntensity> &v){return sum + v.size();}
		);


	const QVector<int> trancheSizes = calculateTranchSize(pointCountAllScans, threadCount);
	e = ErrorUtils::isEqual(trancheSizes.size(), threadCount); ree;

	scanNumberMzIntensityTranched->resize(threadCount);
	int pointCounter = 0;
	int currentIndex = 0;
	int currentTargetTrancheSize = trancheSizes[currentIndex];
	for (QVector<ScanNumberMzIntensity> &v: m_precursorMzKeyVsScanNumberMzIntensity) {
		for (ScanNumberMzIntensity &scanNumberMzIntensity : v) {
			if (pointCounter++ == currentTargetTrancheSize) {
				currentIndex++;
				currentTargetTrancheSize = trancheSizes[currentIndex];
				pointCounter = 0;
			}

			e = ErrorUtils::isTrue(currentIndex < threadCount); ree;
			(*scanNumberMzIntensityTranched)[currentIndex].push_back(&scanNumberMzIntensity);
		}
	}

	// int counter = 0;
	// for (const QVector<ScanNumberMzIntensity*> &v : scanNumberMzIntensityTranched) {
	// 	qDebug()
	// 	<< qPrintable(S_GLOBAL_TIMER.elapsed())
	// 	<< counter++
	// 	<< "Size" << v.size()
	// 	<< "Precursor Mz range"
	// 	<< v.front()->scanInfoPntr->precursorTargetMz
	// 	<< "-"
	// 	<< v.back()->scanInfoPntr->precursorTargetMz
	// 	;
		// for (auto a : v) {
		// 	std::cout << a->scanInfoPntr->precursorTargetMz << " " << a->mzVal << std::endl;
		// }
	// }

	ERR_RETURN
}
