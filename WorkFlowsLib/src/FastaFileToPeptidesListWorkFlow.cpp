//
// Created by anichols on 4/4/23.
//

#include "FastaFileToPeptidesListWorkFlow.h"

#include "BiophysicalCalcs.h"
#include "DIAMzTargetsReader.h"
#include "ErrorUtils.h"
#include "FastaReader.h"
#include "MathUtils.h"
#include "ProteinDigestomatic.h"
#include "TandemFragmentPredictotron.h"

#include <QtConcurrent/QtConcurrent>
#include <QFuture>

#include <iostream>
#include <random>


Err FastaFileToPeptidesListWorkFlow::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}


namespace {

    Err writePeptideListToParquet(
            const QString &fastaFilePath,
            const QVector<PeptideSequence> &peptideSequences,
            QString *outputFilePath
            ) {

        ERR_INIT

        *outputFilePath = fastaFilePath + S_GLOBAL_SETTINGS.DOT_PEPLIB;

        const QVector<QSharedPointer<ParquetReaderInputBase>> peptideSequencesPointers
                = ParquetReaderInputBase::convertInputStructToSharedPointers(peptideSequences);

        ParquetReader parquetReader;
        e = parquetReader.writeDataToParquet(
                *outputFilePath,
                peptideSequencesPointers
        ); ree;

        qDebug() << "Peptides generate:" << peptideSequences.size();
        qDebug() << "Peptide Library written to:" << *outputFilePath;

        ERR_RETURN
    }


}//namespace
Err FastaFileToPeptidesListWorkFlow::exec(
        const QString &fastaFilePath,
        const QString &targetMzCollisionCSV,
        QString *outputFilePath
        ) {

    ERR_INIT

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;

    const QMap<ProteinId, FastaEntry> fastaEntries = fastaReader.fastaEntries();

    QVector<PeptideSequence> peptideSequences;
    e = digestFastaEntries(
            fastaEntries,
            &peptideSequences
            ); ree;

    if (m_params.addDecoys) {
        e = addDecoys(
                666,
                &peptideSequences
                ); ree;
    }

    e = writeLibraryBuilderCSV(
            peptideSequences,
            targetMzCollisionCSV
            ); ree;

    e = writePeptideListToParquet(
            fastaFilePath,
            peptideSequences,
            outputFilePath
            ); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<PeptideSequence>> digestLogic(
            const FastaEntry &fe,
            const PythiaParameters &params
            ) {

        ERR_INIT

        ProteinDigestomatic proteinDigestomatic(params);
        QVector<PeptideSequence> peptideSequences;

        e = proteinDigestomatic.digestProtein(
                fe.fastaSequence,
                &peptideSequences
                );

        return {e, peptideSequences};
    }

}//namespace
Err FastaFileToPeptidesListWorkFlow::digestFastaEntries(
        const QMap<ProteinId, FastaEntry> &fastaEntries,
        QVector<PeptideSequence> *peptideSequences
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fastaEntries); ree;

    const auto digestLogicBinder = std::bind(
            digestLogic,
            std::placeholders::_1,
            m_params
    );

    QFuture<QPair<Err, QVector<PeptideSequence>>> futures = QtConcurrent::mapped(
            fastaEntries,
            digestLogicBinder
            );
    futures.waitForFinished();

    QHash<PeptideString , bool> entered;
    for (const QPair<Err, QVector<PeptideSequence>> &result : futures) {
        e = result.first; ree;

        for (const PeptideSequence &ps : result.second) {

            if (entered.value(ps.sequence)) {
                continue;
            }

            entered[ps.sequence] = true;
            peptideSequences->push_back(ps);
        }
    }

    qDebug() << "Targets size" << peptideSequences->size();
    ERR_RETURN
}


namespace {

    const QString FAILED_SHUFFLE = QStringLiteral("Failed Shuffle");

    QString shufflePeptide(const QString& peptideSeq) {

        std::mt19937 rng(666); //TODO make this settable.

        const int minSizePeptideShuffle = 5;
        const int peptideLen = peptideSeq.size();
        if (peptideLen < minSizePeptideShuffle) {
            return peptideSeq;
        }

        QStringList strList = peptideSeq.split("", QString::SkipEmptyParts);

        std::reverse(strList.rbegin() + 1, strList.rend() - 1);

        std::shuffle(
                strList.begin() + 1,
                strList.begin() + peptideLen - 1,
                rng
        );

        return strList.join("");
    }

    PeptideSequence parallelAddDecoysLogic(const PeptideSequence &ps) {

        PeptideSequence newPepSeq = ps;
        newPepSeq.isDecoy = true;

        int retries = 0;
        const int maxRetryCount = 10;
        while (newPepSeq.sequence == ps.sequence) {

            newPepSeq.sequence = shufflePeptide(ps.sequence);

            if (retries++ > maxRetryCount) {
                newPepSeq.sequence = FAILED_SHUFFLE;
                return newPepSeq;
            }
        }

        return newPepSeq;
    }

}//namespace
Err FastaFileToPeptidesListWorkFlow::addDecoys(
        int seed,
        QVector<PeptideSequence> *peptideSequences
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*peptideSequences); ree;
    const int peptideSequencesSize = peptideSequences->size();

#define RUN_PARALLEL_DECOYS
#ifdef RUN_PARALLEL_DECOYS

    QFuture<PeptideSequence> futures = QtConcurrent::mapped(
            *peptideSequences,
            parallelAddDecoysLogic
    );
    futures.waitForFinished();

    for (const PeptideSequence &pepSeq : futures) {

        if (pepSeq.sequence == FAILED_SHUFFLE) {
            continue;
        }

        peptideSequences->append(pepSeq);
    }

    e = ErrorUtils::isTrue(peptideSequencesSize != peptideSequences->size()); ree;
    qDebug() << "Targets + Decoys size" << peptideSequences->size();

#else
    QVector<PeptideSequence> newPepSeqs;

    for (const PeptideSequence &ps : *peptideSequences) {
        const PeptideSequence newPepSeq = parallelAddDecoysLogic(ps);

        if (newPepSeq.sequence == FAILED_SHUFFLE) {
            continue;
        }

        newPepSeqs.push_back(newPepSeq);
    }

    peptideSequences->append(newPepSeqs);
#endif

    ERR_RETURN
}


namespace {

    Err buildtargetMzVsCollisionEnergyLookUpMap(
            const QString &targetMzCollisionCSV,
            QMap<double, double> *lookupMap,
            double *maxAlloableMz
            ) {

        ERR_INIT

        DIAMzTargetsReader diaMzTargetsReader;
        QVector<DIAMzTargetsReaderRow> diaMzTargetRows;

        e = diaMzTargetsReader.read(
                targetMzCollisionCSV,
                &diaMzTargetRows
        ); ree;

        e = ErrorUtils::isNotEmpty(diaMzTargetRows); ree;

        const auto sortTargetMzAscLogic
            = [](const DIAMzTargetsReaderRow &l, const DIAMzTargetsReaderRow &r){
            return l.targetMz < r.targetMz;
        };

        std::sort(diaMzTargetRows.begin(), diaMzTargetRows.end(), sortTargetMzAscLogic);

        const double buffer = 1.0;
        for (const DIAMzTargetsReaderRow &r : diaMzTargetRows) {
            lookupMap->insert(r.targetMz, r.collisionEnergy);
            *maxAlloableMz = r.targetMz + r.isoWindowHi + buffer;
        }

        ERR_RETURN
    }

}//namespace
Err FastaFileToPeptidesListWorkFlow::writeLibraryBuilderCSV(
        const QVector<PeptideSequence> &peptideSequences,
        const QString &targetMzCollisionCSV
) {

    ERR_INIT

    double maxAllowableMz = -1.0;
    QMap<double, double> targetMzVsCollisionEnergy;
    e = buildtargetMzVsCollisionEnergyLookUpMap(
            targetMzCollisionCSV,
            &targetMzVsCollisionEnergy,
            &maxAllowableMz
            ); ree;

    const QVector<double> ceLookUpVecKey = targetMzVsCollisionEnergy.keys().toVector();
    const QVector<double> ceLookUpVecVals = targetMzVsCollisionEnergy.values().toVector();

    QString outputFilePath = targetMzCollisionCSV;
    outputFilePath = outputFilePath.replace(
            S_GLOBAL_SETTINGS.DOT_CSV,
            S_GLOBAL_SETTINGS.DOT_LIB + S_GLOBAL_SETTINGS.DOT_CSV
    );

    QVector<PeptidePredictionInput> rowsToWrite;

    for (int chrg = m_params.chargeStateMin; chrg <= m_params.chargeStateMax; ++chrg) {

        for (const PeptideSequence &ps : peptideSequences) {

            const double mz = BiophysicalCalcs::calculateThomsonFromMass(ps.mass, chrg);
            const int ceLookUpVecValsIndex = MathUtils::closest(ceLookUpVecKey, mz);
            const double collisionEnergy = ceLookUpVecVals.at(ceLookUpVecValsIndex);

            if (mz > maxAllowableMz) {
                continue;
            }

            PeptidePredictionInput row;
            row.peptideSequence = ps.sequence;
            row.charge = chrg;
            row.collisionEnergy = collisionEnergy;

            rowsToWrite.push_back(row);
        }
    }

    QVector<QSharedPointer<CSVReaderInputBase>> sharedPtrs
            = CSVReaderInputBase::convertInputStructToSharedPointers(rowsToWrite); ree;

    CSVReader reader;
    e = reader.writeDataToCSV(
            outputFilePath,
            sharedPtrs
            ); ree;

    qDebug() << "Library build csv written to:" << outputFilePath;
    ERR_RETURN
}
