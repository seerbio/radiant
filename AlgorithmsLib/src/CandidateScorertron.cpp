//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "ErrorUtils.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include "EigenUtils.h"

#include <Eigen/Dense>


CandidateScorertron::CandidateScorertron(int topNMS2Ions) : m_topNMS2Ions(topNMS2Ions) {}


namespace {

    Err buildMzHashedVsIonPresence(const QMap<MzHashed, XICPoints> &xicPointMap) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(xicPointMap); ree;

        for (auto it = xicPointMap.begin(); it != xicPointMap.end(); it++) {
            qDebug() << it.key();
            qDebug() << it.value().scanNumbersVsIntensityVals;
        }

        ERR_RETURN
    }

    int getMaxRowCount(
            const QVector<MS2Ion> &ms2Ions,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence
    ) {

        int maxRowCount = 0;
        for (const MS2Ion &ms2Ion : ms2Ions) {

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            maxRowCount = std::max(maxRowCount, ionPresenceVec.size());
        }

        return maxRowCount;
    }

    Eigen::MatrixX<double> buildSummingMatrix(
            const QVector<MS2Ion> &ms2Ions,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            int topNMs2Ions
    ) {

        const int cols = topNMs2Ions;
        const int rows = getMaxRowCount(
                ms2Ions,
                mzHashedVsIonPresence
        );

        if (rows == 0) {
            return {};
        }

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            if (ionPresenceVec.isEmpty()) {
                continue;
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(ionPresenceVec);
            mat.col(i) = eVec;
        }

        return mat;
    }
}
Err CandidateScorertron::calculateScores(
        const QMap<MzHashed, XICPoints> &xicPointMap,
        const QVector<MS2Ion> &ms2Ions
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xicPointMap); ree;

    e = buildMzHashedVsIonPresence(xicPointMap); ree;

//    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
//            ms2Ions,
//            m_mzHashedVsIonPresence,
//            m_topNMS2Ions
//    );


    ERR_RETURN
}


