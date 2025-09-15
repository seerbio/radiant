//
// Created by anichols on 9/26/23.
//

#include "FDRCLassifierNeuralNet.h"

#include "BiophysicalCalcs.h"
#include "CandidateScores.h"
#include "CandidateScoresDDA.h"
#include "EigenUtils.h"
#include "GlobalSettings.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"

#include <random>


FDRCLassifierNeuralNet::FDRCLassifierNeuralNet()
: m_epochs(-1)
, m_batchSize(-1)
, m_learningRate(-1.0)
, m_isInit(false)
, m_threadCount(8)
, m_baggingSize(6)
, m_nodesFraction(0.5)
{}

FDRCLassifierNeuralNet::~FDRCLassifierNeuralNet() {
    for (CandidateClassifier *cc : m_candidateClassifiers) {
        delete cc;
    }
}

Err FDRCLassifierNeuralNet::init(
        int epochs,
        int baggingSize,
        int batchSize,
        double learningRate,
        double nodesFraction,
        int threadCount
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(epochs > 0); ree;
    e = ErrorUtils::isTrue(baggingSize >= 1); ree;
    e = ErrorUtils::isTrue(batchSize > 0); ree;
    e = ErrorUtils::isTrue(threadCount > 0); ree;
    e = ErrorUtils::isTrue(learningRate > 0 & learningRate < 1 & !MathUtils::tZero(learningRate)); ree;

    m_epochs = epochs;
    m_baggingSize = baggingSize;
    m_batchSize = batchSize;
    m_learningRate = learningRate;
    m_threadCount = threadCount;
    m_nodesFraction = nodesFraction;

    m_isInit = true;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::exec(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed,
        int verbosity,
        QVector<float> *meanPredictions
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(xData); ree;

    e = trainClassifier(
            xData,
            yData,
            seed,
            verbosity
            ); ree;

    e = predictBaggedClassifiers(xData, meanPredictions); ree;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::trainClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed,
        int verbosity
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData); ree;
    e = ErrorUtils::isEqual(xData.size(), yData.size());

    e = trainBaggedNeuralNets(
            xData,
            yData,
            seed,
            verbosity
            ); ree;

    ERR_RETURN
}

namespace {

    struct CandidateClassifierParallelInput {
        CandidateClassifier *candidateClassifier = nullptr;
        QVector<QVector<float>> xData;
        QVector<float> yData;
        int epochs = -1;
        int batchSize = -1;
        double learningRate = -1.0;
        double nodesFraction = -1.0;
        int bag = -1;
    };

    Err trainingLogic(
        const CandidateClassifierParallelInput &input,
        int verbosity
        ) {

        ERR_INIT

        QMutex mutex;
        mutex.lock();

        e = ErrorUtils::isNotEmpty(input.xData); ree;
        e = ErrorUtils::isEqual(input.xData.size(), input.yData.size()); ree;

        bool trainingCompletedNoErrors = input.candidateClassifier->trainCandidateClassifier(
                input.xData,
                input.yData,
                input.epochs,
                input.batchSize,
                input.learningRate,
                input.bag,
                input.nodesFraction,
                verbosity
        );
        e = ErrorUtils::isTrue(trainingCompletedNoErrors); ree;

        mutex.unlock();

        ERR_RETURN
    }

}//namespace
Err FDRCLassifierNeuralNet::trainBaggedNeuralNets(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed,
        int verbosity
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData);
    e = ErrorUtils::isEqual(xData.size(), yData.size()); ree;

    QElapsedTimer et;
    et.start();

    QVector<CandidateClassifierParallelInput> parallelInputs;
    for (int bag = 0; bag < m_baggingSize; bag++) {

        auto *candidateClassifier = new CandidateClassifier();
        m_candidateClassifiers.push_back(candidateClassifier);

        CandidateClassifierParallelInput ccpi;
        ccpi.candidateClassifier = candidateClassifier;
        ccpi.xData = xData;
        ccpi.yData = yData;
        ccpi.epochs = m_epochs;
        ccpi.batchSize = m_batchSize;
        ccpi.learningRate = m_learningRate;
        ccpi.nodesFraction = m_nodesFraction;
        ccpi.bag = bag + S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST + seed;

        parallelInputs.push_back(ccpi);
    }

// #define PARALLEL_TRAIN_NN
#ifdef PARALLEL_TRAIN_NN
    const auto trainingLogicBinder = std::bind(
        trainingLogic,
        std::placeholders::_1,
        verbosity
        );

    QFuture<Err> futures = QtConcurrent::mapped(
            parallelInputs,
            trainingLogicBinder
            );
    futures.waitForFinished();

    for (const Err &result : futures) {
        e = result; ree;
    }
#else
    for (const CandidateClassifierParallelInput &ccpi : parallelInputs) {
        e = trainingLogic(ccpi, verbosity); ree;
    }
#endif

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Neural nets trained in:"
            << et.elapsed()
            << "mSec";
    ERR_RETURN
}

Err FDRCLassifierNeuralNet::predictBaggedClassifiers(
        const QVector<QVector<float>> &allDataVecs,
        QVector<float> *meanPredictions
        ) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_candidateClassifiers); ree;
    e = ErrorUtils::isNotEmpty(allDataVecs); ree;

    Eigen::VectorX<double> predictionsSum(allDataVecs.size());ree
    predictionsSum.setZero();

    for (const CandidateClassifier *cc : m_candidateClassifiers) {

        QVector<float> predictions;
        const bool predictionOK = cc->predict(
                allDataVecs,
                &predictions
                );
        e = ErrorUtils::isTrue(predictionOK); ree;

        QVector<double> predictionsDouble(predictions.begin(), predictions.end());

        const Eigen::VectorX<double> tranchePred = EigenUtils::convertQVectorToEigenVector(predictionsDouble);
        predictionsSum += tranchePred;
    }

    predictionsSum /= static_cast<double>(m_baggingSize);

    const QVector<double> vD = EigenUtils::convertEigenVectorToQVector(predictionsSum);

    *meanPredictions = QVector<float>(vD.begin(), vD.end());

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        const QVector<CandidateScores> &candidateScores,
        double qValueThreshold,
        int *targetCountBelowFDRThreshold
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    const auto countLogic = [qValueThreshold](const CandidateScores &cs){
        return cs.qValue < qValueThreshold;
    };

    *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(candidateScores.begin(), candidateScores.end(), countLogic));

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        QVector<CandidateScores*> &candidateScores,
        double qValueThreshold,
        int *targetCountBelowFDRThreshold
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    const auto countLogic = [qValueThreshold](const CandidateScores *cs){
        return cs->qValue < qValueThreshold;
    };

    *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(candidateScores.begin(), candidateScores.end(), countLogic));

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
	QVector<CandidateScoresDDA*> &candidateScores,
	double qValueThreshold,
	int *targetCountBelowFDRThreshold
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(candidateScores); ree;
	e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

	const auto countLogic = [qValueThreshold](const CandidateScoresDDA *cs){
		return cs->qValue < qValueThreshold;
	};

	*targetCountBelowFDRThreshold
			= static_cast<int>(std::count_if(candidateScores.begin(), candidateScores.end(), countLogic));

	ERR_RETURN

}

Err FDRCLassifierNeuralNet::outputFDRResults(
        const QVector<CandidateScores> &candidateScores,
        bool verbose,
        QMap<int, int> *fdrVsCount
        ) {

    ERR_INIT

    const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01};
    for (double fdrThresh : fdrFractions) {
        int foundAtThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                candidateScores,
                fdrThresh,
                &foundAtThreshold
        ); ree;
        const int fdrPercent = static_cast<int>(fdrThresh * 100);
        fdrVsCount->insert(fdrPercent, foundAtThreshold);

        if (!verbose) {
            continue;
        }

        QString builder;
        for (auto it = fdrVsCount->begin(); it != fdrVsCount->end(); ++it) {
            const QString &k = QString::number(it.key());
            builder += k + "%: " + QString::number(it.value()) + " | ";
        }

        if (verbose > 0) {
            qDebug() << "PSMs found:" << builder;
        }
    }

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::outputFDRResults(
        QVector<CandidateScores*> &candidateScores,
        int verbose,
        QMap<int, int> *fdrVsCount
        ) {

    ERR_INIT

    const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01};
    for (double fdrThresh : fdrFractions) {
        int foundAtThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                candidateScores,
                fdrThresh,
                &foundAtThreshold
        ); ree;
        const int fdrPercent = static_cast<int>(fdrThresh * 100);
        fdrVsCount->insert(fdrPercent, foundAtThreshold);
    }

    QString outputString;
    e = outPutFDRCounts(*fdrVsCount, &outputString); ree;

    if (verbose > 0) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMs found:" << qPrintable(outputString);
    }

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::outputFDRResults(
	QVector<CandidateScoresDDA*> &candidateScores, int verbose,
	QMap<int, int> *fdrVsCount
	) {

	ERR_INIT

   const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01};
	for (double fdrThresh : fdrFractions) {
		int foundAtThreshold;
		e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
				candidateScores,
				fdrThresh,
				&foundAtThreshold
		); ree;
		const int fdrPercent = static_cast<int>(fdrThresh * 100);
		fdrVsCount->insert(fdrPercent, foundAtThreshold);
	}

	QString outputString;
	e = outPutFDRCounts(*fdrVsCount, &outputString); ree;

	if (verbose > 0) {
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMs found:" << qPrintable(outputString);
	}

	ERR_RETURN
}

Err FDRCLassifierNeuralNet::outPutFDRCounts(
    const QMap<int, int> &fdrVsCount,
    QString *outputString
    ) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(fdrVsCount); ree;

    QStringList builder;
    for (auto it = fdrVsCount.begin(); it != fdrVsCount.end(); ++it) {
        const QString &k = QString::number(it.key());
        builder.push_back(k + "%: " + QString::number(it.value()));
    }

    *outputString = builder.join(" | ");

    ERR_RETURN
}
