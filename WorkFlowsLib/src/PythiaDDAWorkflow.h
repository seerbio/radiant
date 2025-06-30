//
// Created by andrewnichols on 6/29/25.
//

#ifndef PYTHIADDAWORKFLOW_H
#define PYTHIADDAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

class MsReaderMzMLLazyLoad;
class MsScanInfo;

using namespace Error;


class WORKFLOWSLIB_EXPORTS PythiaDDAWorkflow {

public:

	PythiaDDAWorkflow() = default;
	~PythiaDDAWorkflow() = default;

	Err init();

	Err processFile(const QString &msDataFilePath);


private:

	Err processChunk(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		size_t msScanInfosMinIndex,
		size_t msScanInfosMaxIndex,
		MsReaderMzMLLazyLoad *msReaderMzMlLazyLoad
		);


private:

};



#endif //PYTHIADDAWORKFLOW_H
