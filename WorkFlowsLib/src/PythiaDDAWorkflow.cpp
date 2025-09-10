//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DeIsotopotron.h"
#include "EigenUtils.h"
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
	m_parameters.ms2ExtractionWidthPPM = 7;
	m_parameters.threadCount = 64; //TODO use higher threadcount to load balans
	m_parameters.print();

	e = FragLibReader::getFragLibReaderRows(
		libraryFilePath,
		&m_fragLibReaderRows
		); ree;

	e = m_tdcpManager.init(
		m_parameters,
		&m_fragLibReaderRows
		); ree;

	m_tdcpManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyCandidatePairsPntrs);

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

	e = m_msFraggertron.init(m_parameters, &msReaderPtr); ree;
	e = buildMsCalibratomatic(); ree;


	ERR_RETURN
}

Err PythiaDDAWorkflow::buildMsCalibratomatic() {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairsPntrs); ree;

	QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsCalibration;
	for (int i = 0; i < m_targetDecoyCandidatePairsPntrs.size(); i++) {
		constexpr int skipCount = 1;
		if (i % skipCount != 0 ) {
			continue;
		}
		targetDecoyCandidatePairsCalibration.push_back(m_targetDecoyCandidatePairsPntrs[i]);
	}

	e = m_msFraggertron.performFragging(targetDecoyCandidatePairsCalibration); ree;


	ERR_RETURN
}
