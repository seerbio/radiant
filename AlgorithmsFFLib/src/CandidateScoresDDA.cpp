//
// Created by ubuntu on 9/9/25.
//

#include "CandidateScoresDDA.h"

void CandidateScoresDDA::initFeaturesArray() {
	featuresArray = QVector<float>(FeaturesDDASize, 0.0f);
	featuresArray.reserve(FeaturesDDASize);
}

QVector<float> CandidateScoresDDA::selectFeaturesArrayFeatures(const QVector<FeaturesDDA> &enumFeatures) {

	QVector<float> selectedFeatures(enumFeatures.size());

	for (int i = 0; i < enumFeatures.size(); i++) {
		selectedFeatures[i] = featuresArray[enumFeatures.at(i)];
	}

	return selectedFeatures;
}

QVector<float> CandidateScoresDDA::selectFeaturesArrayFeatures(
		const QVector<float> &featureArray,
		const QVector<FeaturesDDA> &enumFeatures
		) {
	QVector<float> selectedFeatures(enumFeatures.size());

	for (int i = 0; i < enumFeatures.size(); i++) {
		selectedFeatures[i] = featureArray[enumFeatures.at(i)];
	}

	return selectedFeatures;

}

void CandidateScoresDDA::printFeatures(const QVector<FeaturesDDA> &featuresToPrint) {

	QVector<float> vec;
	for (const FeaturesDDA &f : featuresToPrint) {
		vec.push_back(featuresArray[f]);
	}
	qDebug() << vec;
}