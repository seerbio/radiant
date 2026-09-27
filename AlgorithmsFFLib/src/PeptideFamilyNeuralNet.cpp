#include "PeptideFamilyNeuralNet.h"

#include "CandidateClassifier.h"
#include "DiscriminantScoretron.h"
#include "EigenUtils.h"

#include <QCryptographicHash>
#include <QFuture>
#include <QRegularExpression>
#include <QSet>
#include <QThreadPool>
#include <QtConcurrent/QtConcurrentRun>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <numeric>
#include <random>

namespace {

bool ionMobilityOnly(Features feature) {
    switch (feature) {
    case Ms1IntensityFoundApex100IM:
    case IonMobilityDelta:
    case IonMobilityDeltaAbs:
    case IonMobilityPdAbs:
    case Ms2IonMobilityWeightedDelta:
    case Ms2IonMobilityWeightedDeltaAbs:
    case Ms2IonMobilityApexDeltaAbsMean:
    case Ms2IonMobilityApexDeltaAbsStDev:
    case Ms2IonMobilityMatchedIonFraction:
    case Ms2IonMobilityFwhmMean:
    case Ms2IonMobilityFwhmStDev:
    case Ms2IonMobilityRtCosineMean:
    case Ms2IonMobilityRtCosineStDev:
    case Ms2IonMobilityRtApexAgreementFraction:
        return true;
    default:
        return false;
    }
}

bool rawIntensity(Features feature) {
    return feature == TotalIntensityRaw || feature == Ms1IntensityFound100
        || feature == Ms1IntensityFound45 || feature == Ms1IntensityFoundPreMono
        || feature == Ms1IntensityFoundIso1 || feature == Ms1IntensityFoundIso2
        || feature == Ms1IntensityFoundApex100IM;
}

} // namespace

quint32 PeptideFamilyNeuralNet::familyHash(QString originSequence) {
    originSequence.remove(QRegularExpression("\\([^)]*\\)"));
    originSequence.remove('_');
    originSequence.replace('I', 'L');
    return QCryptographicHash::hash(originSequence.toUtf8(), QCryptographicHash::Sha256)
        .toHex().left(8).toUInt(nullptr, 16);
}

QVector<Features> PeptideFamilyNeuralNet::nonTimsFeatures() {
    auto features = DiscriminantScoretron::featuresNeuralNetwork();
    features.erase(std::remove_if(features.begin(), features.end(), ionMobilityOnly), features.end());
    return features;
}

Error::Err PeptideFamilyNeuralNet::score(
    const QVector<QVector<float>> &rawFeatures,
    const QVector<quint32> &decoyLabels,
    const QVector<quint32> &familyHashes,
    const QVector<Features> &features,
    const Settings &settings,
    Predictions *predictions) {
    const int count = rawFeatures.size();
    if (predictions == nullptr || count < 2 || count != decoyLabels.size()
        || count != familyHashes.size() || features.isEmpty()
        || settings.folds < 2 || settings.folds > 64 || settings.folds > count
        || settings.networks < 1 || settings.networks > 64
        || settings.epochs < 1 || settings.threads < 1 || settings.threads > 64
        || settings.seed < 0 || settings.seed > (std::numeric_limits<int>::max() - settings.networks) / 2
        || !std::isfinite(settings.learningRate) || settings.learningRate <= 0 || settings.learningRate >= 1
        || !std::isfinite(settings.nodesFraction) || settings.nodesFraction <= 0 || settings.nodesFraction > 4
        || features.size() * settings.nodesFraction < 1.0
        || !std::isfinite(settings.focalLossGamma) || settings.focalLossGamma < 0
        || qint64(count) * settings.networks > std::numeric_limits<int>::max()) {
        return Error::eValueError;
    }
    QSet<int> featureIds;
    for (Features feature : features) {
        if (feature < 0 || feature >= FeaturesSize || featureIds.contains(int(feature)))
            return Error::eValueError;
        featureIds.insert(int(feature));
    }
    if (!settings.normalizationGroups.isEmpty()) {
        qint64 rows = 0;
        for (int group : settings.normalizationGroups) {
            if (group < 1) return Error::eValueError;
            rows += group;
        }
        if (rows != count) return Error::eValueError;
    }

    QVector<QVector<float>> values;
    values.reserve(count);
    for (int row = 0; row < count; ++row) {
        if (rawFeatures[row].size() != FeaturesSize || decoyLabels[row] > 1)
            return Error::eValueError;
        auto selected = CandidateScores::selectFeaturesArrayFeatures(rawFeatures[row], features);
        for (int column = 0; column < selected.size(); ++column) {
            if (!std::isfinite(selected[column])) return Error::eValueError;
            if (settings.logIntensities && rawIntensity(features[column]))
                selected[column] = std::log1p(std::max(0.0f, selected[column]));
        }
        values.push_back(std::move(selected));
    }
    auto matrix = EigenUtils::convertQVectorsToEigenMatrix(values);
    if (settings.normalizationGroups.isEmpty()) {
        EigenUtils::minMaxScaleMatrix(&matrix);
    } else {
        int offset = 0;
        for (int rows : settings.normalizationGroups) {
            Eigen::MatrixX<float> group = matrix.middleRows(offset, rows);
            EigenUtils::minMaxScaleMatrix(&group);
            matrix.middleRows(offset, rows) = group;
            offset += rows;
        }
    }
    values = EigenUtils::convertEigenMatrixToQVectors(matrix);
    for (const auto &row : values)
        for (float value : row)
            if (!std::isfinite(value)) return Error::eValueError;

    QVector<int> shuffled(count);
    std::iota(shuffled.begin(), shuffled.end(), 0);
    std::mt19937 generator(settings.seed);
    for (int pass = 0; pass < 3; ++pass)
        std::shuffle(shuffled.begin(), shuffled.end(), generator);
    QVector<QVector<int>> members(settings.folds);
    for (int row : shuffled)
        members[familyHashes[row] % settings.folds].push_back(row);
    for (const auto &fold : members)
        if (fold.isEmpty() || fold.size() == count) return Error::eValueError;
    const auto &normalized = values;
    const auto &groups = members;

    using FoldResult = QPair<Error::Err, QVector<float>>;
    const auto fit = [&](int heldOut) -> FoldResult {
        try {
            QVector<QVector<float>> training, inference;
            QVector<float> labels;
            training.reserve(count - groups[heldOut].size());
            labels.reserve(count - groups[heldOut].size());
            inference.reserve(groups[heldOut].size());
            for (int fold = 0; fold < settings.folds; ++fold) {
                for (int row : groups[fold]) {
                    if (fold == heldOut) {
                        inference.push_back(normalized[row]);
                    } else {
                        training.push_back(normalized[row]);
                        labels.push_back(float(decoyLabels[row]));
                    }
                }
            }
            QVector<float> result(inference.size() * settings.networks);
            for (int network = 0; network < settings.networks; ++network) {
                CandidateClassifier classifier;
                if (!classifier.trainCandidateClassifier(
                        training, labels, settings.epochs,
                        std::min(500, std::max(2, count / 100)),
                        settings.learningRate, 2 * settings.seed + network,
                        settings.nodesFraction, settings.focalLossGamma, 0,
                        settings.shuffleEachEpoch)) {
                    return {Error::eError, {}};
                }
                QVector<float> output;
                if (!classifier.predict(inference, &output) || output.size() != inference.size())
                    return {Error::eError, {}};
                for (int row = 0; row < output.size(); ++row) {
                    if (!std::isfinite(output[row]) || output[row] < 0 || output[row] > 1)
                        return {Error::eValueError, {}};
                    result[row * settings.networks + network] = output[row];
                }
            }
            return {Error::eNoError, result};
        } catch (const std::exception &error) {
            qWarning() << "Family neural-net training failed:" << error.what();
            return {Error::eError, {}};
        }
    };
    // A private pool avoids changing the caller's global concurrency settings.
    QThreadPool pool;
    pool.setMaxThreadCount(std::min(settings.threads, settings.folds));
    QVector<QFuture<FoldResult>> futures;
    for (int fold = 0; fold < settings.folds; ++fold)
        futures.push_back(QtConcurrent::run(&pool, std::function<FoldResult()>(
            [&, fold] { return fit(fold); })));
    pool.waitForDone();

    Predictions result;
    result.perNetwork.resize(count * settings.networks);
    result.meanDecoyProbability.resize(count);
    result.heldOutFold.resize(count);
    for (int fold = 0; fold < settings.folds; ++fold) {
        const auto fitted = futures[fold].result();
        if (fitted.first != Error::eNoError) return fitted.first;
        for (int index = 0; index < members[fold].size(); ++index) {
            const int row = members[fold][index];
            double total = 0.0;
            for (int network = 0; network < settings.networks; ++network) {
                const float probability = fitted.second[index * settings.networks + network];
                result.perNetwork[row * settings.networks + network] = probability;
                total += probability;
            }
            result.meanDecoyProbability[row] = float(total / settings.networks);
            result.heldOutFold[row] = quint32(fold);
        }
    }
    *predictions = std::move(result);
    return Error::eNoError;
}
