//
// Created by andrewnichols on 5/23/25.
//

#ifndef MSFRAGGATRON_H
#define MSFRAGGATRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"

using namespace Error;

class TargetDecoyCandidatePair;

class ALGORITHMSFFLIB_EXPORTS MsFraggertron {

public:
	MsFraggertron();
	~MsFraggertron();

	Err init(QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> *mzHashedVsTDCPPntrs);

	Err findTargetDecoyCandidates(
		double mz,
		double ppmTolerance,
		QVector<TargetDecoyCandidatePair*> *targetDecoyCandidates
		);


private:

	Q_DISABLE_COPY(MsFraggertron) class Private;
	const QScopedPointer<Private> d_ptr;

};




#endif //MSFRAGGATRON_H
