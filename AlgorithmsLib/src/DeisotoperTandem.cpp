//
// Created by Drucifer on 4/27/2022.
//

#include "DeisotoperTandem.h"

#include <QtConcurrent/QtConcurrent>
#include "EigenSparseUtils.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "ParallelUtils.h"

#include <vector>


namespace {

    Eigen::SparseVector<double> buildSubtractionVector(
            const ScanPoints &scanPoints,
            int precision,
            double ppmTol
            ) {

        const auto sortLogic = [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();};
        const double maxMz = std::max_element(scanPoints.begin(), scanPoints.end(), sortLogic)->x();

        const int buffer = 10;
        const int hashedMaxXVal = MathUtils::hashDecimal(maxMz + buffer, precision);

        Eigen::SparseVector<double> vec(hashedMaxXVal);
        vec.setZero();

        const double subtractionFactor = 3.0;

        for (const ScanPoint &p : scanPoints) {

            const double mzTol = MathUtils::calculatePPM(p.x(), ppmTol);
            const int mzTolHashed = std::max(MathUtils::hashDecimal(mzTol, precision), 1);

            const int mzHashed = MathUtils::hashDecimal(p.x(), precision);
            const int startMzHashed = mzHashed - mzTolHashed;
            const int endMzHashed = mzHashed + mzTolHashed + 1;

            const double yValBuffered = p.y() * subtractionFactor;

            for (int insertMzHashed = startMzHashed; insertMzHashed < endMzHashed; ++insertMzHashed) {

                if (insertMzHashed < 0 || insertMzHashed >= hashedMaxXVal) {
                    continue;
                }

                vec.coeffRef(insertMzHashed) = yValBuffered;
            }

        }

        return vec;
    }

}//namespace
ScanPoints DeisotoperTandem::deisotopeTandemScan(
        const ScanPoints &tandemScan,
        double ppmTolerance
        ) {

    const auto sortLogic = [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();};
    const double maxMz = std::max_element(tandemScan.begin(), tandemScan.end(), sortLogic)->x();

    Eigen::SparseVector<double> vec
        = EigenSparseUtils::vectorizeFromScan(tandemScan, maxMz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

    const Eigen::SparseVector<double> subtractionVecOG = buildSubtractionVector(
            tandemScan,
            S_GLOBAL_SETTINGS.HASHING_PRECISION,
            ppmTolerance
            );

    int currentCharge = S_GLOBAL_SETTINGS.MAX_CHARGE_TANDEM_DEISOTOPING;
    while (currentCharge > 0) {

        const double rollDistance = S_GLOBAL_SETTINGS.ISO_DIFF / currentCharge ;
        const int rollDistanceHashed = MathUtils::hashDecimal(rollDistance, S_GLOBAL_SETTINGS.HASHING_PRECISION);

        const Eigen::SparseVector<double> subtractionVecRolled
            = EigenSparseUtils::roll(subtractionVecOG, rollDistanceHashed);

        vec -= subtractionVecRolled;
        --currentCharge;
    }

    const double thresholdValue = 0.0;
    EigenSparseUtils::removeElementsBelowThreshold(thresholdValue, &vec);

    ScanPoints returnVec;
    returnVec.reserve(static_cast<int>(vec.size()));

    for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it) {

        const auto mz
            = MathUtils::unHashDecimal<double>(static_cast<int>(it.index()), S_GLOBAL_SETTINGS.HASHING_PRECISION);

        const double intensity = it.value();

        returnVec.push_back({mz, intensity});
    }

    return returnVec;
}


namespace {

    QPair<ScanNumber, ScanPoints> parallelDeisotopeLogic(
            const QPair<ScanNumber, ScanPoints> &pr,
            double ppmTolerance
            ) {

        return {pr.first,
                DeisotoperTandem::deisotopeTandemScan(
                        pr.second,
                        ppmTolerance
                        )
        };
    }


}//namespace
Err DeisotoperTandem::deisotopeTandemScansParallel(
        const QMap<ScanNumber, ScanPoints> &tandemScans,
        double ppmTolerance,
        QMap<ScanNumber, ScanPoints> *tandemScansDeisotoped
        ) {

    ERR_INIT

    tandemScansDeisotoped->clear();
    e = ErrorUtils::isNotEmpty(tandemScans); ree;

    const QVector<QPair<ScanNumber, ScanPoints>> scanPointsAsPair
            = ParallelUtils::convertMapToVectorPairs<ScanNumber, ScanPoints>(tandemScans);

    const auto deisotoperLogicBinder = std::bind(
            parallelDeisotopeLogic,
            std::placeholders::_1,
            ppmTolerance
    );

    QFuture<QPair<ScanNumber, ScanPoints>> futures = QtConcurrent::mapped(
            scanPointsAsPair,
            deisotoperLogicBinder
            );
    futures.waitForFinished();

    for (const QPair<ScanNumber, ScanPoints> &pr : futures) {
        tandemScansDeisotoped->insert(pr.first, pr.second);
    }

    ERR_RETURN
}

