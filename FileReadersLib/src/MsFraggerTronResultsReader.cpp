//
// Created by anichols on 2/22/23.
//

#include "MsFraggerTronResultsReader.h"

#include "ErrorUtils.h"
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

namespace {

     Err decodeQStringToDoubleMz(
             const QString &stringOfDoubles,
             QVector<double> *vec
                     ) {

         ERR_INIT

        const QStringList stringOfDoublesSplit
            = stringOfDoubles.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

         for (const QString &mzStr : stringOfDoublesSplit) {

            double mz;
            e = ErrorUtils::toDouble(mzStr, &mz); ree;
            vec->push_back(mz);
        }

        ERR_RETURN
    }

}//namespace
Err MsFraggerTronResultsReader::readCsv(
        const QString &firstPassCsvFilePath,
        QVector<RowToWrite> *rowsToWrite
        ) {

    ERR_INIT

    QFile file(firstPassCsvFilePath);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << file.errorString();
        rrr(eFileError);
    }

    const QString expectedHeader = "ScanNumber,PeptideId,Occurrence,IntensityTotal,MeanMzPPM,"
                                   "Sequence,SequenceWithMods,PreviousResidue,PostResidue,"
                                   "Mass,IsDecoy,ScanIonMz,TheoFragIonMz\n";

    QString header;
    while (!file.atEnd()) {

        if (header.isEmpty()) {
            header = file.readLine();

            e = ErrorUtils::isEqualString(header, expectedHeader); ree;
            continue;
        }

        const QString &row = file.readLine().replace(S_GLOBAL_SETTINGS.NEWLINE, "");
        const QStringList rowSplit = row.split(S_GLOBAL_SETTINGS.COMMA);

        RowToWrite rtw;

        e = ErrorUtils::toInt(rowSplit.at(0), &rtw.scanNumber); ree;
        e = ErrorUtils::toInt(rowSplit.at(1), &rtw.peptideId); ree;
        e = ErrorUtils::toInt(rowSplit.at(2), &rtw.occurrence); ree;
        e = ErrorUtils::toInt(rowSplit.at(3), &rtw.intensityTotal); ree;
        e = ErrorUtils::toDouble(rowSplit.at(4), &rtw.meanMzPPM); ree;
        rtw.sequence = rowSplit.at(5);
        rtw.sequenceWithMods = rowSplit.at(6);
        rtw.previousResidue = rowSplit.at(7).at(0);
        rtw.postResidue = rowSplit.at(8).at(0);
        e = ErrorUtils::toDouble(rowSplit.at(9), &rtw.mass); ree;

        int boolConversion;
        e = ErrorUtils::toInt(rowSplit.at(10), &boolConversion); ree;
        rtw.isDecoy = boolConversion;

        e = decodeQStringToDoubleMz(rowSplit.at(11), &rtw.scanIonMZs); ree;
        e = decodeQStringToDoubleMz(rowSplit.at(12), &rtw.theoFragIonMZs); ree;

        rowsToWrite->push_back(rtw);
    }

    ERR_RETURN
}
