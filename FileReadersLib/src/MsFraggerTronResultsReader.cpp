//
// Created by anichols on 2/22/23.
//

#include "MsFraggerTronResultsReader.h"

#include "SqlUtils.h"

#include <QFile>

#include <iostream>

namespace {

    void sortRowsToWriteOccurrenceDescending(QVector<RowToWrite> *rowsToWrite) {

        const auto sortLogic = [](const RowToWrite &l, const RowToWrite &r){
            return l.occurrence < r.occurrence;
        };

        std::sort(rowsToWrite->rbegin(), rowsToWrite->rend(), sortLogic);
    }
}
void MsFraggerTronResultsReader::writeToCsv(
        const QString &outputFilePath,
        QVector<RowToWrite> *rowsToWrite
        ) {

    sortRowsToWriteOccurrenceDescending(rowsToWrite);

    const QStringList rowHeader = {
            QStringLiteral("ScanNumber"),
            QStringLiteral("PeptideId"),
            QStringLiteral("Occurrence"),
            QStringLiteral("IntensityTotal"),
            QStringLiteral("MeanMzPPM"),
            QStringLiteral("Sequence"),
            QStringLiteral("SequenceWithMods"),
            QStringLiteral("PreviousResidue"),
            QStringLiteral("PostResidue"),
            QStringLiteral("Mass"),
            QStringLiteral("IsDecoy"),
            QStringLiteral("ScanIonMz"),
            QStringLiteral("TheoFragIonMz")
    };

    int counter = 0;

    QFile file(outputFilePath);
    if (file.open(QIODevice::ReadWrite)) {

        QTextStream stream(&file);
        stream << rowHeader.join(S_GLOBAL_SETTINGS.COMMA) + S_GLOBAL_SETTINGS.NEWLINE;

        for (const RowToWrite &rtw : *rowsToWrite) {

            const QStringList row = {
                    QString::number(rtw.scanNumber),
                    QString::number(rtw.peptideId),
                    QString::number(rtw.occurrence),
                    QString::number(rtw.intensityTotal),
                    QString::number(rtw.meanMzPPM),
                    rtw.sequence,
                    rtw.sequenceWithMods,
                    rtw.previousResidue,
                    rtw.postResidue,
                    QString::number(rtw.mass),
                    QString::number(rtw.isDecoy)
            };
            stream << row.join(S_GLOBAL_SETTINGS.COMMA) + S_GLOBAL_SETTINGS.COMMA;

            for (double mz : rtw.scanIonMZs) {
                stream << QString::number(mz) << S_GLOBAL_SETTINGS.SEPARATOR;
            }

            stream << S_GLOBAL_SETTINGS.COMMA;

            for (double mz : rtw.theoFragIonMZs) {
                stream << QString::number(mz) << S_GLOBAL_SETTINGS.SEPARATOR;
            }

            stream << S_GLOBAL_SETTINGS.NEWLINE;

            if (counter++ == 1000) { //TODO make sure to delete this.
                break;
            }
        }
    }

    file.close();
    qDebug() << "PSM results written to:" << outputFilePath;
}
