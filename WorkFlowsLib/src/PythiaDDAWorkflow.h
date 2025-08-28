//
// Created by andrewnichols on 6/29/25.
//

#ifndef PYTHIADDAWORKFLOW_H
#define PYTHIADDAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

class MsReaderPointerAcc;
class MsScanInfo;

using namespace Error;


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

	Err performFragging(const QString &msDataFilePath);


private:

	PythiaParameters m_parameters;
	QVector<FragLibReaderRow> m_fragLibReaderRows;
	TargetDecoyCandidatePairManager m_tdcpManager;
	QVector<QVector<MsScanInfo*>> m_msScanInfosTranched;
	QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsPntrs;
	QMap<ScanNumber, MsScanInfo*> m_scanNumberVsMsScanInfoMS2;

};



#endif //PYTHIADDAWORKFLOW_H
