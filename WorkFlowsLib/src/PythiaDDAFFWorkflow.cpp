//
// Created by andrewnichols on 5/22/25.
//

#include "PythiaDDAFFWorkflow.h"

#include "FragLibReader.h"
#include "MsReaderPointerAcc.h"


PythiaDDAFFWorkflow::PythiaDDAFFWorkflow() {}
PythiaDDAFFWorkflow::~PythiaDDAFFWorkflow() {}


Err PythiaDDAFFWorkflow::init(
    const PythiaParameters &pythiaParameters,
    const QString &fragLibUri,
    const QString &fastaUri,
    const QString &outputFolderPath
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    ParallelUtils::printSystemDetails();

    m_fastaUri = fastaUri;
    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;
    m_outputFolderPath = outputFolderPath;
    m_pythiaParameters.print();

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading library";

    e = FragLibReader::getFragLibReaderRows(
            m_fragLibUri,
            &m_fragLibReaderRows
            ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished reading library";

    e = m_targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            &m_fragLibReaderRows
            ); ree;

    e = buildMzHashedVsTargetDecoyCandidatePairs(); ree;
	e = m_msFraggertron.init(&m_mzHashedVsTargetDecoyCandidatePairs); ree;

    ERR_RETURN
}

namespace {

	struct DDAJobChunk {
		PythiaParameters pythiaParameters;
		int scanNumberMin = -1;
		int scanNumberMax = -1;
		MsReaderPointerAcc *msReaderPointerAcc = nullptr;
		MsFraggertron *msFraggertron = nullptr;
	};

	QVector<DDAJobChunk> buildDdaJobChunks(
		const PythiaParameters &pythiaParameters,
		int ms2ScanCount,
		int chunkSize,
		MsFraggertron *msFraggertron,
		MsReaderPointerAcc *msReaderPointerAcc
	) {
		QVector<DDAJobChunk> jobChunks;

		for (int i = 0; i < ms2ScanCount; i += chunkSize) {
			DDAJobChunk chunk;
			chunk.pythiaParameters = pythiaParameters;
			chunk.scanNumberMin = i + 1; //scanNumber does not start at 0
			chunk.scanNumberMax = std::min(chunk.scanNumberMin + chunkSize - 1, ms2ScanCount);
			chunk.msFraggertron = msFraggertron;
			chunk.msReaderPointerAcc = msReaderPointerAcc;
			jobChunks.push_back(chunk);
		}

		return jobChunks;
	}

	Err processDDAJobChunkLogic(const DDAJobChunk &chunk) {

		ERR_INIT

		e = ErrorUtils::isTrue(chunk.msReaderPointerAcc->isInit()); ree;
		e = ErrorUtils::isTrue(chunk.scanNumberMin < chunk.scanNumberMax); ree;

		QVector<QPair<MsScanInfo*, ScanPoint*>> scanNumberVsScanPointPntrs;
		e = chunk.msReaderPointerAcc->ptr->buildMsScanInfoPntrVsScanPointPntrsMs2(
			chunk.scanNumberMin,
			chunk.scanNumberMax,
			chunk.pythiaParameters.threadCount,
			&scanNumberVsScanPointPntrs
			); ree;

		if (chunk.pythiaParameters.verbosity > -1) { //TODO change verbosity threshold from -1 -> 0
			qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "processDDAJobChunkLogic" << chunk.scanNumberMin;
		}


		ERR_RETURN
	}

}//namespace
Err PythiaDDAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;

    msReaderPointerAcc.setUseLazyLoading(false);

    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    msReaderPointerAcc.ptr->printSize();

	constexpr int chunkSize = 1000; //TODO find a better way to set this based on scancounts
	const int lastScanNumber = msReaderPointerAcc.ptr->lastScanNumber();

	const QVector<DDAJobChunk> ddaJobChunks = buildDdaJobChunks(
		m_pythiaParameters,
		lastScanNumber,
		chunkSize,
		&m_msFraggertron,
		&msReaderPointerAcc
		);

	QFuture<Err> future = QtConcurrent::mapped(
		ddaJobChunks,
		processDDAJobChunkLogic
		);
	future.waitForFinished();



    ERR_RETURN
}

namespace {

    QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> buildMzHashedVsTargetDecoyCandidatePairsLogic(
        const QVector<TargetDecoyCandidatePair*> &tdcps
        ) {

        QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> mzHashedVsTargetDecoyCandidatePairs;

        for (TargetDecoyCandidatePair* tdcp : tdcps) {
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsDecoy()) {
                const auto mz = static_cast<double>(ms2Ion.mz);
                const int mzHashed = MathUtils::hashDecimal(mz, S_GLOBAL_SETTINGS.HASHING_PRECISION_DDA);
                mzHashedVsTargetDecoyCandidatePairs[mzHashed].push_back(tdcp);
            }
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsTarget()) {
                const auto mz = static_cast<double>(ms2Ion.mz);
                const int mzHashed = MathUtils::hashDecimal(mz, S_GLOBAL_SETTINGS.HASHING_PRECISION_DDA);
                mzHashedVsTargetDecoyCandidatePairs[mzHashed].push_back(tdcp);
            }
        }

        return mzHashedVsTargetDecoyCandidatePairs;
    }

}//namespace
Err PythiaDDAFFWorkflow::buildMzHashedVsTargetDecoyCandidatePairs() {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairs;
    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&targetDecoyCandidatePairs); ree;

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        targetDecoyCandidatePairs,
        m_pythiaParameters.threadCount,
        &targetDecoyCandidatePairsTranched
        ); ree;

    QFuture<QHash<MzHashed, QVector<TargetDecoyCandidatePair*>>> futures = QtConcurrent::mapped(
        targetDecoyCandidatePairsTranched,
        buildMzHashedVsTargetDecoyCandidatePairsLogic
        );
    futures.waitForFinished();

    for (const QHash<MzHashed, QVector<TargetDecoyCandidatePair*>> &result : futures) {
        for (auto it = result.begin(); it != result.end(); ++it) {
            const MzHashed mzHashed = it.key();
            const QVector<TargetDecoyCandidatePair*> &tdcps = it.value();
            m_mzHashedVsTargetDecoyCandidatePairs[mzHashed].append(tdcps);
        }
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building" << m_mzHashedVsTargetDecoyCandidatePairs.size() << "hashes";

    ERR_RETURN
}
