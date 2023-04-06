//
// Created by anichols on 4/4/23.
//

#include "FastaFileToPeptidesListWorkFlow.h"

#include "ErrorUtils.h"
#include "FastaReader.h"
#include "ProteinDigestomatic.h"
#include <QtConcurrent/QtConcurrent>
#include <QFuture>


Err FastaFileToPeptidesListWorkFlow::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}


Err FastaFileToPeptidesListWorkFlow::exec(
        const QString &fastaFilePath,
        QString *outputFilePath
        ) {

    ERR_INIT

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;

    const QMap<ProteinId, FastaEntry> fastaEntries = fastaReader.fastaEntries();

    QVector<QSharedPointer<ParquetReaderInputBase>> peptideSequences;
    e = digestFastaEntries(
            fastaEntries,
            &peptideSequences
            ); ree;

    *outputFilePath = fastaFilePath + S_GLOBAL_SETTINGS.DOT_PEPLIB;

    ParquetReader parquetReader;
    e = parquetReader.writeDataToParquet(
            *outputFilePath,
            peptideSequences
            ); ree;

    qDebug() << "Peptides generate:" << peptideSequences.size();
    qDebug() << "Peptide Library written to:" << *outputFilePath;

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
        QVector<QSharedPointer<ParquetReaderInputBase>> *peptideSequences
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fastaEntries); ree;

    const auto degestLogicBinder = std::bind(
            digestLogic,
            std::placeholders::_1,
            m_params
    );

    QFuture<QPair<Err, QVector<PeptideSequence>>> futures = QtConcurrent::mapped(
            fastaEntries,
            degestLogicBinder
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

            QSharedPointer<ParquetReaderInputBase> parquetReaderInputBase(new PeptideSequence(ps));
            peptideSequences->push_back(parquetReaderInputBase);
        }
    }

    ERR_RETURN
}
