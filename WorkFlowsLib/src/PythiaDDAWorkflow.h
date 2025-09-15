//
// Created by andrewnichols on 6/29/25.
//

#ifndef PYTHIADDAWORKFLOW_H
#define PYTHIADDAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFraggertron.h"
#include "MsCalibratomatic.h"
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

	Err buildMsCalibratomatic();


private:


	Err honeIRTAndMassCalibration(
		QVector<CandidateScoresDDA*> *candidateScoresVecScoredPntrs,
		int topNCandidates,
		int topCandidatesMass
	);



private:

	PythiaParameters m_parameters;
	MsFraggertron m_msFraggertron;
	MsCalibratomatic m_msCalibratomatic;
	QVector<FragLibReaderRow> m_fragLibReaderRows;
	TargetDecoyCandidatePairManager m_tdcpManager;
	QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsPntrs;

};

#endif //PYTHIADDAWORKFLOW_H
