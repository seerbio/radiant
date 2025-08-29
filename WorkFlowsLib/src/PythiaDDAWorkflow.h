//
// Created by andrewnichols on 6/29/25.
//

#ifndef PYTHIADDAWORKFLOW_H
#define PYTHIADDAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

class MsReaderPointerAcc;
class MsScanInfo;

struct MsScanPoint {
	MsScanInfo *scanInfoPntr = nullptr;
	float mzVal = -1.0;
	float intensityVal = -1.0;
};

struct MS2IonLibrary {
	TargetDecoyCandidatePair *targeDecoyCandidatePairPntr = nullptr;
	MS2Ion *ms2IonPntr = nullptr;
	bool isDecoy = false;
};

struct ProcessingGroup {
	QPair<float, float> mzPrecursorRangeMinMax;
	QVector<MsScanPoint*> msScanPoints;
	QVector<MS2IonLibrary*> ms2IonsLibrary;
};

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

	Err extractScansParallel(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		);
	Err performFragging();


private:

	PythiaParameters m_parameters;
	QVector<FragLibReaderRow> m_fragLibReaderRows;
	TargetDecoyCandidatePairManager m_tdcpManager;
	QVector<QVector<MsScanInfo*>> m_msScanInfosTranched;
	QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsPntrs;
	QMap<ScanNumber, MsScanInfo*> m_scanNumberVsMsScanInfoMS2;
	QVector<QVector<MsScanPoint>> m_msScanPointsTranched;

};



#endif //PYTHIADDAWORKFLOW_H
