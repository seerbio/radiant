//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DeIsotopotron.h"
#include "MsReaderMzMLLazyLoad.h"
#include "ObjectCSVWriters.h"
#include "ParallelUtils.h"

Err PythiaDDAWorkflow::init() {

	ERR_INIT


	ERR_RETURN
}

Err PythiaDDAWorkflow::processFile(const QString &msDataFilePath) {

	ERR_INIT

	MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
	e = msReaderMzMLLazyLoad.openFile(msDataFilePath); ree;

	constexpr size_t msLevel = 2;

	const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderMzMLLazyLoad.getMsScanInfos(msLevel);
	QVector<MsScanInfo> scanNumbers = msScanInfos.values().toVector();

	QVector<MsScanInfo*> scanInfosPntrs;
	std::transform(
		scanNumbers.begin(),
		scanNumbers.end(),
		std::back_inserter(scanInfosPntrs),
		[](MsScanInfo &scanInfo) {return &scanInfo;}
		);


	constexpr int chunkSize = 1e4;
	for (int chunkIndexMin = 0; chunkIndexMin < scanNumbers.size(); chunkIndexMin += chunkSize) {
		e = processChunk(
			scanInfosPntrs,
			chunkIndexMin,
			std::min(chunkIndexMin + chunkSize, scanNumbers.size()),
			&msReaderMzMLLazyLoad
			); ree;
	}

	ERR_RETURN
}

namespace {

	QPair<Err, QMap<ScanNumber, ScanPoints>> readScanPointsParallel(
		const QVector<MsScanInfo*> &msScanInfosSubset,
		float ppmTol,
		MsReaderMzMLLazyLoad *msReaderMzMlLazyLoad
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(msScanInfosSubset); rree;

		QMap<ScanNumber, ScanPoints> scanPoints;
		e = msReaderMzMlLazyLoad->extractScanPoints(
			msScanInfosSubset,
			&scanPoints
			); rree;

		for (auto it = scanPoints.begin(); it != scanPoints.end(); ++it) {
			const ScanNumber scanNumber = it.key();
			ScanPoints &scanPointsIter = it.value();

			DeIsotopotron::deisotopeTandemScan(ppmTol, &scanPointsIter);
		}

		return {e, scanPoints};
	}

	Err processChunkParallelLogic(
		const QVector<MsScanInfo*> &msScanInfosSubset,
		MsReaderMzMLLazyLoad *msReaderMzMlLazyLoad
		) {

		ERR_INIT

		QVector<QVector<MsScanInfo*>> msScanInfosSubsetTranched;
		e = ParallelUtils::trancheVectorForParallelization(
			msScanInfosSubset,
			ParallelUtils::numberOfAvailableSystemProcessors(),
			&msScanInfosSubsetTranched
			); ree;

// #define DDA_PARALLEL
#ifdef DDA_PARALLEL

		const auto binderLogic = std::bind(
			readScanPointsParallel,
			std::placeholders::_1,
			10.0, //TODO replace this w/ extractionPPM
			msReaderMzMlLazyLoad
			);

		QFuture<QPair<Err, QMap<ScanNumber, ScanPoints>>> futures = QtConcurrent::mapped(
			msScanInfosSubsetTranched,
			binderLogic
			);
		futures.waitForFinished();

		for (const QPair<Err, QMap<ScanNumber, ScanPoints>> &result : futures) {
			e = result.first; ree;
		}
#else

		QPair<Err, QMap<ScanNumber, ScanPoints>> result = readScanPointsParallel(
			msScanInfosSubset,
			10.0f,
			msReaderMzMlLazyLoad
			);
		e = result.first; ree;


#endif



		ERR_RETURN
	}

}//namespace
Err PythiaDDAWorkflow::processChunk(
	const QVector<MsScanInfo*> &scanInfosPntrs,
	size_t msScanInfosMinIndex,
	size_t msScanInfosMaxIndex,
	MsReaderMzMLLazyLoad *msReaderMzMlLazyLoad
	) {

	ERR_INIT

	const QVector<MsScanInfo*> msScanInfosSubset(
		scanInfosPntrs.begin() + msScanInfosMinIndex,
		scanInfosPntrs.begin() + msScanInfosMaxIndex
		);


	e = processChunkParallelLogic(
		msScanInfosSubset,
		msReaderMzMlLazyLoad
		); ree;



	ERR_RETURN
}







