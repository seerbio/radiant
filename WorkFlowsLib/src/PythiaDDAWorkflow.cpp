//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DeIsotopotron.h"
#include "MsReaderMzMLLazyLoad.h"
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

	QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs;
	m_tdcpManager.getTargetDecoyCandidatePairPointers(&targetDecoyCandidatePairsPntrs);

	qDebug() << targetDecoyCandidatePairsPntrs.size() << "SDJFL" << m_fragLibReaderRows.size();


	ERR_RETURN
}

Err PythiaDDAWorkflow::processFile(const QString &msDataFilePath) {

	ERR_INIT

	MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
	e = msReaderMzMLLazyLoad.openFile(msDataFilePath); ree;

	// constexpr size_t msLevel = 2;
	//
	// const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderMzMLLazyLoad.getMsScanInfos(msLevel);
	// QVector<MsScanInfo> scanNumbers = msScanInfos.values().toVector();
	//
	// QVector<MsScanInfo*> scanInfosPntrs;
	// std::transform(
	// 	scanNumbers.begin(),
	// 	scanNumbers.end(),
	// 	std::back_inserter(scanInfosPntrs),
	// 	[](MsScanInfo &scanInfo) {return &scanInfo;}
	// 	);
	//
	//
	// constexpr int chunkSize = 1e4;
	// for (int chunkIndexMin = 0; chunkIndexMin < scanNumbers.size(); chunkIndexMin += chunkSize) {
	// 	e = processChunk(
	// 		scanInfosPntrs,
	// 		chunkIndexMin,
	// 		std::min(chunkIndexMin + chunkSize, scanNumbers.size()),
	// 		&msReaderMzMLLazyLoad
	// 		); ree;
	// }

	ERR_RETURN
}




