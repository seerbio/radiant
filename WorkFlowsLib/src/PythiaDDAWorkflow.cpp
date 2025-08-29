//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DeIsotopotron.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
#include "ParallelUtils.h"

Err PythiaDDAWorkflow::init(
	const PythiaParameters &parameters,
	const QString& libraryFilePath
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(parameters.isValid()); ree;
	e = ErrorUtils::fileExists(libraryFilePath); ree;

	m_parameters = parameters;
	m_parameters.print();

	e = FragLibReader::getFragLibReaderRows(
		libraryFilePath,
		&m_fragLibReaderRows
		); ree;

	e = m_tdcpManager.init(
		m_parameters,
		&m_fragLibReaderRows
		); ree;

	int idCounter = 0;
	m_tdcpManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyCandidatePairsPntrs);
	for (TargetDecoyCandidatePair *tdcp : m_targetDecoyCandidatePairsPntrs) {
		tdcp->id = idCounter;
		idCounter += 2;
	}

	qDebug()
	<< qPrintable(S_GLOBAL_TIMER.elapsed())
	<< m_targetDecoyCandidatePairsPntrs.size()
	<< "targets from library of" << m_fragLibReaderRows.size()
	<< "after filter.";

	ERR_RETURN
}

Err PythiaDDAWorkflow::processFile(const QString &msDataFilePath) {

	ERR_INIT

	MsReaderPointerAcc msReaderPtr;
	e = msReaderPtr.openFile(msDataFilePath); ree;

	constexpr int msLevel = 2;
	m_scanNumberVsMsScanInfoMS2 = msReaderPtr.ptr->getMsScanInfos(msLevel); ree;
	e = ErrorUtils::isNotEmpty(m_scanNumberVsMsScanInfoMS2); ree;

	e = extractScansParallel(
		m_scanNumberVsMsScanInfoMS2.values().toVector(),
		&msReaderPtr
		); ree;

	// e = performFragging(); ree;

	ERR_RETURN
}

namespace {

	QPair<Err, QVector<MsScanPoint>> extractScansParallelLogic(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		) {
		ERR_INIT

		QVector<MsScan> msScans;
		e = msReaderPtr->ptr->extractScanPoints(
			scanInfosPntrs,
			&msScans
			); rree;

		QVector<MsScanPoint> msScanPoints;
		for (const MsScan &scan : msScans) {
			for (int i = 0; i < scan.msScanInfoPntr->pointCount; i++) {
				MsScanPoint msScanPoint;
				msScanPoint.scanInfoPntr = scan.msScanInfoPntr;
				msScanPoint.mzVal = scan.mzVals[i];
				msScanPoint.intensityVal = scan.intensityVals[i];
				msScanPoints.push_back(msScanPoint);
			}

		}

		constexpr float mzMin = 300; //TODO get values from MsReader for these vals
		constexpr float mzMax = 1600; //TODO get values from MsReader for these vals
		const auto terminatorLogic = [mzMin, mzMax](const MsScanPoint &mssp) {
			return mssp.mzVal < mzMin || mssp.mzVal > mzMax;
		};
		const auto terminator = std::remove_if(msScanPoints.begin(), msScanPoints.end(), terminatorLogic);
		msScanPoints.erase(terminator, msScanPoints.end());

		std::sort(
			msScanPoints.begin(),
			msScanPoints.end(),
			[](const MsScanPoint &l, const MsScanPoint &r) {
				if (MathUtils::tSame(l.scanInfoPntr->precursorTargetMz, r.scanInfoPntr->precursorTargetMz)) {
					return l.mzVal < r.mzVal;
				}
				return l.scanInfoPntr->precursorTargetMz < r.scanInfoPntr->precursorTargetMz;
			});

		return {e, msScanPoints};
	}

}//namespace
Err PythiaDDAWorkflow::extractScansParallel(
	const QVector<MsScanInfo*> &scanInfosPntrs,
	MsReaderPointerAcc *msReaderPtr
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;

	QVector<QVector<MsScanInfo*>> scanInfosPntrsTranched;
	e = ParallelUtils::trancheVectorForParallelization(
		scanInfosPntrs,
		m_parameters.threadCount * 4,
		&scanInfosPntrsTranched
		); ree;

	const auto binderLogic = std::bind(
		extractScansParallelLogic,
		std::placeholders::_1,
		msReaderPtr
		);

	QFuture<QPair<Err, QVector<MsScanPoint>>> future = QtConcurrent::mapped(
		scanInfosPntrsTranched,
		binderLogic
		);
	future.waitForFinished();

	for (const QPair<Err, QVector<MsScanPoint>> &res : future) {
		e = res.first; ree;
		m_msScanPointsTranched.push_back(res.second);
	}

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsPoints extracted";

	ERR_RETURN
}

namespace {

	std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> buildMs2Ions(
		TargetDecoyCandidatePair *tdcp
		) {

		ERR_INIT

		constexpr float mzMin = 300; //TODO get values from MsReader for these vals
		constexpr float mzMax = 1600; //TODO get values from MsReader for these vals
		const QVector<MS2Ion> &ms2IonsTarget = tdcp->ms2IonsTarget(mzMin, mzMax);
		const QVector<MS2Ion> &ms2IonsDecoy = tdcp->ms2IonsDecoy(ms2IonsTarget);

		std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> res(
			tdcp, ms2IonsTarget, ms2IonsDecoy
			);

		return res;
	}

	QPair<Err, int> performFraggingLogic(
		const QVector<MsScan*> &msScansPntrs,
		const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions,
		const PythiaParameters &parameters
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(ms2Ions); rree;
		e = ErrorUtils::isNotEmpty(msScansPntrs); rree;

		Ms2IonFraggertronManager fragger;
		e = fragger.init(ms2Ions); rree;

		for (const MsScan *msScan : msScansPntrs) {

			const MzVals &mzVals = msScan->mzVals ;
			const IntensityVals &intensityVals = msScan->intensityVals ;
			const ScanNumber scanNumber = msScan->msScanInfoPntr->scanNumber;
			const int pointCount = msScan->msScanInfoPntr->pointCount;

			float cutoffIntensity = -1.0;
			constexpr int topNScanPoints = 500;
			if (pointCount > topNScanPoints) {
				QVector<float> intensityValsCopy(intensityVals.data(), intensityVals.data() + pointCount);
				std::sort(intensityValsCopy.rbegin(), intensityValsCopy.rend());
				cutoffIntensity = intensityValsCopy[topNScanPoints - 1];
			}
			for (int i = 0; i < pointCount; ++i) {

				const float intensity = msScan->intensityVals[i];
				if (intensity < cutoffIntensity) {
					continue;
				}

				constexpr float iRTMin = -10000;
				constexpr float iRTMax = 10000;

				const float mzVal = msScan->mzVals[i];
				const float mzTol = MathUtils::calculatePPM(mzVal, static_cast<float>(parameters.ms2ExtractionWidthPPM));
				const float mzMin = mzVal - mzTol;
				const float mzMax = mzVal + mzTol;

				QVector<MS2Frag*> tdPeptideFrags;
				e = fragger.extractMs2Points(
					mzMin,
					mzMax,
					iRTMin,
					iRTMax,
					&tdPeptideFrags
					); rree;
			}
		}

		return {e, -1};
	}


}//namespace
Err PythiaDDAWorkflow::performFragging() {

	ERR_INIT

	QElapsedTimer et;
	et.start();

	QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
	constexpr int libTrancheSize = 4;
	e = ParallelUtils::trancheVectorForParallelization(
		m_targetDecoyCandidatePairsPntrs,
		libTrancheSize,
		&targetDecoyCandidatePairsPntrsTranched
		); ree;



	// QVector<float> mzValsSample;
	// for (const QVector<MsScan*> &msScans : msScansPntrsTranched) {
	// 	constexpr int skipCount = 50;
	// 	int counter = 0;
	// 	for (const MsScan *msScan : msScans) {
	// 		if (counter++ % skipCount != 0) {
	// 			continue;
	// 		}
	// 		const QVector<float> mzVals(
	// 			msScan->mzVals.data(),
	// 			msScan->mzVals.data() + msScan->msScanInfoPntr->pointCount
	// 			);
	//
	// 		mzValsSample.append(mzVals);
	// 	}
	// }
	// std::sort(mzValsSample.begin(), mzValsSample.end());
	// qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Sample mz distribution size" << mzValsSample.size();
	//
	// const int increment = mzValsSample.size() / m_parameters.threadCount;
	// for (int i = 0; i < m_parameters.threadCount; i++) {
	// 	qDebug() << i << mzValsSample[std::min(i * increment, mzValsSample.size() - 1)] << "upper limit";
	// }
	//
	//


	// for (const QVector<MsScan*> &msScans : msScansPntrsTranched) {
	//
	// 	QElapsedTimer et2;
	// 	et2.start();
	// 	const int totalPointCount = std::accumulate(
	// 		msScans.begin(),
	// 		msScans.end(),
	// 		0,
	// 		[](int sum, const MsScan *s){return sum + s->msScanInfoPntr->pointCount;}
	// 		);
	//
	// 	QVector<MsScanPoint> msScanPoints;
	// 	msScanPoints.reserve(totalPointCount);
	// 	for (MsScan *msScan : msScans) {
	// 		for (int i = 0; i < msScan->msScanInfoPntr->pointCount; ++i) {
	// 			MsScanPoint msScanPoint;
	// 			msScanPoint.mzVal = msScan->mzVals[i];
	// 			msScanPoint.scanInfoPntr = msScan->msScanInfoPntr;
	// 			msScanPoints.push_back(msScanPoint);
	// 		}
	// 	}
	//
	// 	std::sort(
	// 		msScanPoints.begin(),
	// 		msScanPoints.end(),
	// 		[](const MsScanPoint &l, const MsScanPoint &r){return l.mzVal < r.mzVal;}
	// 		);
	//
	// 	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << et2.elapsed() << "SDFLJDS";
	// }


// 	for (const QVector<TargetDecoyCandidatePair*> &tdcps : targetDecoyCandidatePairsPntrsTranched) {
//
// 		QFuture<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> futureLib = QtConcurrent::mapped(
// 			tdcps,
// 			buildMs2Ions
// 			);
// 		futureLib.waitForFinished();
//
// 		QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> ms2Ions;
// 		ms2Ions.reserve(tdcps.size());
// 		for (const std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> &r : futureLib) {
// 			ms2Ions.push_back(r);
// 		}
//
// // #define FRAG_PARALLEL
// // #ifdef FRAG_PARALLEL
// // 		const auto binderLogic = std::bind(
// // 			performFraggingLogic,
// // 			std::placeholders::_1,
// // 			ms2Ions,
// // 			m_parameters
// // 			);
// //
// // 		QFuture<QPair<Err, int>> futureScans = QtConcurrent::mapped(
// // 			msScansPntrsTranched,
// // 			binderLogic
// // 			);
// // 		futureScans.waitForFinished();
// // #else
// // 		for (const QVector<MsScan*> &msScans : msScansPntrsTranched) {
// // 			const QPair<Err, int> res = performFraggingLogic(msScans, ms2Ions, m_parameters);
// // 			e = res.first; ree;
// // 		}
// // #endif
// //

	// }


	ERR_RETURN
}




