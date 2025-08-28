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

	QVector<MsScan> msScans;
	e = MsReaderMzMLLazyLoad::extractScanPoints(
		msDataFilePath,
		m_scanNumberVsMsScanInfoMS2.values().toVector(),
		&msScans
		); ree;

	e = performFragging(msDataFilePath); ree;

	ERR_RETURN
}

namespace {

	using MS2IonsTarget = QVector<MS2Ion>;
	using MS2IonsDecoy = QVector<MS2Ion>;

	std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> buildMs2Ions(
		TargetDecoyCandidatePair *tdcp
		) {

		ERR_INIT

		constexpr float mzMin = 200;
		constexpr float mzMax = 2000;
		const QVector<MS2Ion> &ms2IonsTarget = tdcp->ms2IonsTarget(mzMin, mzMax);
		const QVector<MS2Ion> &ms2IonsDecoy = tdcp->ms2IonsDecoy(ms2IonsTarget);

		std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> res(
			tdcp, ms2IonsTarget, ms2IonsDecoy
			);

		return res;
	}

	QPair<Err, int> performFraggingLogic(
		const QVector<MsScanInfo*> &msScanInfosTranche,
		const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairsPntrs,
		const QString &msDataFilePath
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(msScanInfosTranche); rree;
		e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairsPntrs); rree;
		e = ErrorUtils::fileExists(msDataFilePath); rree;

		QVector<MsScan> msScans;
		e = MsReaderMzMLLazyLoad::extractScanPoints(
			msDataFilePath,
			msScanInfosTranche,
			&msScans
			); rree;

		e = ErrorUtils::isEqual(msScans.size(), msScanInfosTranche.size()); rree;

		constexpr int tdcpTrancheSize = 64;
		QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
		e = ParallelUtils::trancheVectorForParallelization(
			targetDecoyCandidatePairsPntrs,
			tdcpTrancheSize,
			&targetDecoyCandidatePairsPntrsTranched
			); rree;

		for (const QVector<TargetDecoyCandidatePair*> &tdcpPntrs : targetDecoyCandidatePairsPntrsTranched) {
			Ms2IonFraggertronManager fragger;
			e = fragger.init(tdcpPntrs); rree;
		}


		return {e, -1};
	}


}//namespace
Err PythiaDDAWorkflow::performFragging(const QString &msDataFilePath) {

	ERR_INIT

	QElapsedTimer et;
	et.start();

	QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
	constexpr int libTrancheSize = 10;
	e = ParallelUtils::trancheVectorForParallelization(
		m_targetDecoyCandidatePairsPntrs,
		libTrancheSize,
		&targetDecoyCandidatePairsPntrsTranched
		); ree;

	for (const QVector<TargetDecoyCandidatePair*> &tdcps : targetDecoyCandidatePairsPntrsTranched) {

		QFuture<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> future
																= QtConcurrent::mapped(tdcps, buildMs2Ions);
		future.waitForFinished();

		for (const QVector<MsScanInfo*> &msScanInfosTranch : m_msScanInfosTranched) {
			// 	const QPair<Err, int> result = performFraggingLogic(
			// 		msScanInfosTranch,
			// 		m_targetDecoyCandidatePairsPntrs,
			// 		msDataFilePath
			// 		); ree;
			// 	qDebug() << et.restart() << "SDLKFDSJL";
			}
	}


	ERR_RETURN
}




