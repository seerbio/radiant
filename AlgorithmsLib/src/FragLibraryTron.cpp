//
// Created by anichols on 2/19/23.
//

#include "FragLibraryTron.h"

#include "ErrorUtils.h"
#include "MathUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>

namespace {

    void sortFragLibIons(QVector<FragLibIon> &fragLibIons) {

        const double fudgeFactor = 0.000001;
        const auto sortLogic = [fudgeFactor](const FragLibIon &l, const FragLibIon &r){

            const int lmz = static_cast<int>(std::round(l.mzFrag / fudgeFactor));
            const int rmz = static_cast<int>(std::round(r.mzFrag / fudgeFactor));

            if (lmz == rmz) {
                return l.peptideMass < r.peptideMass;
            }

            return lmz < rmz;
        };

        std::sort(fragLibIons.begin(), fragLibIons.end(), sortLogic);
    }

    QVector<FragLibIon> sortFragLibIonsInChunks(const QVector<FragLibIon> &fragLibIons) {

        using NominalMzFrag = int;

        QMap<NominalMzFrag, QVector<FragLibIon>> sortingChunks;
        for (const FragLibIon fli : fragLibIons) {
            const auto nominalMzFrag = static_cast<NominalMzFrag>(std::round(fli.mzFrag));
            sortingChunks[nominalMzFrag].push_back(fli);
        }

        // TODO properly parallelize this.
//        QFuture<void> futures = QtConcurrent::map(
//                sortingChunks,
//                &sortFragLibIons
//        );

        QVector<FragLibIon> sortedFragLibIons;
        for (auto it = sortingChunks.begin(); it != sortingChunks.end(); it++) {
            QVector<FragLibIon> &fragLibIonsAtNominalMzFrag = it.value();
            sortFragLibIons(fragLibIonsAtNominalMzFrag);

            sortedFragLibIons.append(fragLibIonsAtNominalMzFrag);
        }

        return sortedFragLibIons;
    }

}//namespace
Err FragLibraryTron::writeFragLibIons(
        const QVector<FragLibIon> &fragLibIons,
        const QString &fragLibIonsFilePath
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragLibIons); ree;

    QElapsedTimer et;
    et.start();

    QVector<FragLibIon> _fragLibIons = sortFragLibIonsInChunks(fragLibIons);
    qDebug() << "FragLibIons sorted in" << et.restart() << "mSec";

    QFile writeFile(fragLibIonsFilePath);
    if (!writeFile.open(QIODevice::WriteOnly)) {
        qDebug() << "Could not write to file:" << fragLibIonsFilePath
                 << "Error string:" << writeFile.errorString();
        return Error::eFileError;
    }

    QDataStream out(&writeFile);
    out.setVersion(QDataStream::Qt_5_12);

    out << _fragLibIons;

    qDebug() << "FragLibIons written in" << et.elapsed() << "mSec";
    qDebug() << "FragLibIons library written to:" << fragLibIonsFilePath;
    ERR_RETURN
}

Err FragLibraryTron::readFragLibIons(const QString &fragLibIonsFilePath) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    m_fragLibIons.clear();

    QFile readFile(fragLibIonsFilePath);
    QDataStream in(&readFile);

    in.setVersion(QDataStream::Qt_5_12);

    if (!readFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Could not read the file:" << fragLibIonsFilePath
                 << "Error string:" << readFile.errorString();
        return Error::eFileError;
    }

    in >> m_fragLibIons;

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree;

    qDebug() << "FragLibIons loaded from" << fragLibIonsFilePath
             << "in" <<  et.elapsed() << "mSec";

    ERR_RETURN
}

Err FragLibraryTron::buildFragLibIonsUniqueIndexes() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree;

    for (const FragLibIon &fli : m_fragLibIons) {
        m_fragLibIonsUniquePeptideId[fli.hashedKey()].push_back(fli.peptideId);
    }

    qDebug() << "OG" << m_fragLibIons.size() << "OG UNIQUES" << m_fragLibIonsUniquePeptideId.size();


    ERR_RETURN
}
