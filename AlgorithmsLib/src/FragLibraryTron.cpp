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


FragLibraryTron::FragLibraryTron(const PythiaParameters &pythiaParameters)
: m_pythiaParameters(pythiaParameters){}



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

QVector<FragLibIon> FragLibraryTron::fragLibIons() {
    return m_fragLibIons;
}

Err FragLibraryTron::getFragLibIonTranches(
        const QVector<QPair<MzMin, MzMax>> &tranchLimits,
        QVector<FragLibIonTranche> *fragLibIonTranches
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragLibIons); ree;
    e = ErrorUtils::isNotEmpty(tranchLimits); ree;
    e = ErrorUtils::isFalse(m_pythiaParameters.ms2ExtractionWidthPPM == -1); ree;

    fragLibIonTranches->clear();

    for (const QPair<MzMin, MzMax> &trancheLimit : tranchLimits) {

        QPair<int, int> trancheIndexLimits = getFragLibIonTrancheStartStop(trancheLimit);

    }

    ERR_RETURN
}

namespace {

    int getIndexOfLimit(const QVector<FragLibIon> &vec, double mzLimit) {

        auto it = std::min_element(vec.begin(), vec.end(), [mzLimit] (const FragLibIon &a, const FragLibIon &b) {
            return std::abs(a.mzFrag - mzLimit) <= std::abs(b.mzFrag - mzLimit);
        });

        if(it == vec.end()) {
            return -1;
        }

        return it - vec.begin();
    }

} //namespace
QPair<int, int> FragLibraryTron::getFragLibIonTrancheStartStop(const QPair<MzMin, MzMax> &tranchLimits) {

    QPair<int, int> trancheIndexLimits;

    const double toleranceMultiplier = 5.0;

    const double mzMinTolerance
        = MathUtils::calculatePPM(tranchLimits.first, m_pythiaParameters.ms2ExtractionWidthPPM) * toleranceMultiplier;

    const double mzMaxTolerance
            = MathUtils::calculatePPM(tranchLimits.second, m_pythiaParameters.ms2ExtractionWidthPPM) * toleranceMultiplier;

    const int mzMinIndex = getIndexOfLimit(m_fragLibIons, tranchLimits.first - mzMinTolerance);
    const int mzMaxIndex = getIndexOfLimit(m_fragLibIons, tranchLimits.second + mzMaxTolerance);

    return trancheIndexLimits;
}

