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

struct TallyResult {
	ScanNumber scanNumber = -1;
	Occurrence occurrence = 0;
	QVector<int> ranks;
	// int indexesFoundY = 0;
	// int indexesFoundB = 0;
	// int seqTagLongestY = 0;
	// int seqTagLongestB = 0;
};

using TallyResultTarget = TallyResult;
using TallyResultDecoy = TallyResult;
using TallyResultTuple = std::tuple<TargetDecoyCandidatePair*, QVector<TallyResultTarget>, QVector<TallyResultDecoy>>;
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
	QVector<QVector<MsScanPoint*>> m_msScanPointsPntrsTranched;
	QVector<QPair<float, float>> m_mzPrecursorStartVsStopResult;
};



#endif //PYTHIADDAWORKFLOW_H
