//
// Created by andrewnichols on 8/11/25.
//

#include "ApexConnector.h"

#include "EigenUtils.h"

namespace {

	QVector<Eigen::MatrixX<int>> buildApexGroups(const QVector<QVector<int>> &apexes) {

		const int defaultGroupPadding = 2;
		const int nonLastTrancheAdjustment = sgpi.isLastTranche ? 0 : sgpi.skipScanCount;

		for (
				ScanNumberIndex scanIndex = 0;
				scanIndex < sgpi.allScanPointsTranch.size() - (1 + nonLastTrancheAdjustment);
				scanIndex++
		) {

			const ScanNumberIndex stopNextScanIndex
				= scanIndex + sgpi.skipScanCount + defaultGroupPadding;

			QVector<QVector<float>> mzValsGroup;
			QVector<QVector<float>> intensityValsGroup;

			for (ScanNumberIndex nextScanIndex = scanIndex; nextScanIndex < stopNextScanIndex; ++nextScanIndex) {

				if (nextScanIndex >= sgpi.allScanPointsTranch.size()) {
					break;
				}

				ScanPoints* scanPoints = sgpi.allScanPointsTranch.at(nextScanIndex);
				if (scanPoints->isEmpty()) {
					const ScanPoint placeholderPointForEmptyScan(1,0);
					scanPoints->push_back(placeholderPointForEmptyScan);
				}

				QVector<float> mzVals;
				QVector<float> intensityVals;
				output.e = MsReaderBase::splitScanPoints(
						scanPoints,
						&mzVals,
						&intensityVals
				);

				if (output.e != eNoError) {
					return output;
				}

				mzValsGroup.push_back(mzVals);
				intensityValsGroup.push_back(intensityVals);
			}

			output.groupedMzVals.push_back(mzValsGroup);
			output.groupedIntensityVals.push_back(intensityValsGroup);
		}

		return output;

	}

}//namespace
Err ApexConnector::connectApexes(
	const QVector<QVector<int>> &apexes,
	QVector<QVector<int>> *connectedApexes
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(apexes);

	connectedApexes->resize(apexes.size());

	const QVector<Eigen::MatrixX<int>> apexGroups = buildApexGroups(apexes);

	ERR_RETURN
}
