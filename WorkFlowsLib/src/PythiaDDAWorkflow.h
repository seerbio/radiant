//
// Created by andrewnichols on 6/29/25.
//

#ifndef PYTHIADDAWORKFLOW_H
#define PYTHIADDAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderBase.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class MsReaderPointerAcc;

class WORKFLOWSLIB_EXPORTS PythiaDDAWorkflow {

public:

	PythiaDDAWorkflow() = default;
	~PythiaDDAWorkflow() = default;

	Err init(
		const PythiaParameters &parameters,
		const QString& libraryFilePath
		);

	Err processFile(const QString &msDataFilePath);


private:

	Err buildMsScanPointPntrs(MsReaderPointerAcc *msReaderPtr);

	Err extractScansParallel(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		);

	Err performFragging();

	QPair<Err, QVector<TallyResultTuple>> processTargetDecoyCandidatePairsPntrsTranch(
		const QVector<TargetDecoyCandidatePair*> &tdcps,
		const QVector<ProcessingGroup> &processingGroups
		);


private:

	PythiaParameters m_parameters;
	QVector<FragLibReaderRow> m_fragLibReaderRows;
	TargetDecoyCandidatePairManager m_tdcpManager;
	QVector<QVector<MsScanInfo*>> m_msScanInfosTranched;
	QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsPntrs;
	QMap<ScanNumber, MsScanInfo*> m_scanNumberVsMsScanInfoMS2;
	QVector<QVector<MsScanPoint>> m_msScanPointsTranched;
	// QVector<QVector<MsScanPoint*>> m_msScanPointsPntrsTranched;
	QVector<QPair<float, float>> m_mzPrecursorStartVsStopResult;
};



#endif //PYTHIADDAWORKFLOW_H
