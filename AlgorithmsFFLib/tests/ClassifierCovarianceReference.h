#pragma once
// Frozen covariance implementation from 2b6001a8, before the optimization.
// Keep this independent oracle unchanged when editing the production routine.
#include "EigenUtils.h"
#include "ParallelUtils.h"
#include <QtConcurrent/QtConcurrent>

namespace FrozenCovariance {
    struct ParallelLogicInput {
        QVector<QVector<float>*> targets;
        QVector<QVector<float>*> decoys;
        Eigen::VectorX<float> matTargetsSumMean;
        Eigen::VectorX<float> matDecoysSumMean;
    };

    Eigen::MatrixX<float> aMatrixParallelLogic(const ParallelLogicInput &input) {

        const int rows = input.targets.size();
        const int cols = input.targets.front()->size();

        Eigen::MatrixX<double> matA(cols, cols);
        matA.setZero();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                for (int k = j; k < cols; k++) {
                    matA.coeffRef(j, k) += 0.5 *
                                           ((input.decoys[i]->at(j) - input.matDecoysSumMean[j]) * (input.decoys[i]->at(k) - input.matDecoysSumMean[k]) +
                                            (input.targets[i]->at(j) - input.matTargetsSumMean[j]) * (input.targets[i]->at(k) - input.matTargetsSumMean[k]));
                    matA.coeffRef(k, j) = matA.coeffRef(j, k);
                }
            }
        }

        return matA.cast<float>();
    }

    Err buildParallelLogicInput(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        const Eigen::VectorX<float> &matTargetsSumMean,
        const Eigen::VectorX<float> &matDecoysSumMean,
        QVector<ParallelLogicInput> *parallelLogicInputs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targets); ree;
        e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
        e = ErrorUtils::isEqual(targets.front()->size(), decoys.front()->size()); ree;

        const int threadCount = ParallelUtils::numberOfAvailableSystemProcessors();

        QVector<QVector<QVector<float>*>> targetsTranched;
        e = ParallelUtils::trancheVectorForParallelization(
            targets,
            std::min(targets.size(), threadCount),
            &targetsTranched
            ); ree;

        QVector<QVector<QVector<float>*>> decoysTranched;
        e = ParallelUtils::trancheVectorForParallelization(
            decoys,
            std::min(targets.size(), threadCount),
            &decoysTranched
            ); ree;

        e = ErrorUtils::isEqual(targetsTranched.size(), decoysTranched.size()); ree;

        for (int i = 0; i < targetsTranched.size(); i++) {
            ParallelLogicInput parallelLogicInput;
            parallelLogicInput.targets = targetsTranched[i];
            parallelLogicInput.decoys = decoysTranched[i];
            parallelLogicInput.matDecoysSumMean = matDecoysSumMean;
            parallelLogicInput.matTargetsSumMean = matTargetsSumMean;

            parallelLogicInputs->push_back(parallelLogicInput);
        }

        ERR_RETURN
    }

    Err buildDataClassifier2(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        QVector<QVector<float>> *A,
        QVector<float> *b
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targets); ree;
        e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
        e = ErrorUtils::isEqual(targets.front()->size(), decoys.front()->size()); ree;

        const Eigen::MatrixX<float> matTargets = EigenUtils::convertQVectorsToEigenMatrix(targets);
        const Eigen::MatrixX<float> matDecoys = EigenUtils::convertQVectorsToEigenMatrix(decoys);

        const Eigen::VectorX<float> matTargetsSum = matTargets.colwise().sum();
        const Eigen::VectorX<float> matDecoysSum = matDecoys.colwise().sum();

        const Eigen::VectorX<float> matTargetsSumMean = matTargetsSum / matTargets.rows();
        const Eigen::VectorX<float> matDecoysSumMean = matDecoysSum / matTargets.rows();

        const Eigen::MatrixX<float> subtractionMat = matTargets - matDecoys;
        const Eigen::VectorX<float> subtractionMatSum = subtractionMat.colwise().sum();
        const Eigen::VectorX<float> meanMat = subtractionMatSum / matTargets.rows();

        *b = EigenUtils::convertEigenVectorToQVector(meanMat);

        const int rows = targets.size();
        const int cols = targets.front()->size();

        Eigen::MatrixX<float> matA(cols, cols);
        matA.setZero();

        QElapsedTimer et;
        et.start();

#define BUILD_CLASS2_PARALLEL
#ifdef BUILD_CLASS2_PARALLEL

        QVector<ParallelLogicInput> parallelLogicInputs;
        e = buildParallelLogicInput(
            targets,
            decoys,
            matTargetsSumMean,
            matDecoysSumMean,
            &parallelLogicInputs
            ); ree;

        QFuture<Eigen::MatrixX<float>> futures = QtConcurrent::mapped(
            parallelLogicInputs,
            aMatrixParallelLogic
            );
        futures.waitForFinished();

        for (const Eigen::MatrixX<float>& mat : futures) {
            matA += mat;
        }

#else
        ParallelLogicInput input;
        input.targets = targets;
        input.decoys = decoys;
        input.matTargetsSumMean = matTargetsSumMean;
        input.matDecoysSumMean = matDecoysSumMean;
        matA = aMatrixParallelLogic(input);
#endif

        matA /= rows - 1;
        for (int i = 0; i < cols; i++) {
            matA.coeffRef(i, i) += std::numeric_limits<float>::min();
        }

        *A = EigenUtils::convertEigenMatrixToQVectors(matA);

        ERR_RETURN
    }

} // namespace FrozenCovariance
