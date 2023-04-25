//
// Created by anichols on 4/21/23.
//

#include "TandemSpectraDeconvolvotron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "FragLibReader.h"

#include <Eigen/Sparse>


TandemSpectraDeconvolvotron::TandemSpectraDeconvolvotron()
: m_precision(2)
, m_mzMax(1500.0)
, m_iterationsMax(100)
, m_stopTolerance(1e-8)
{}

Err TandemSpectraDeconvolvotron::setParameters(
        int precision,
        double mzMax,
        int iterationsMax,
        double stopTolerance
        ) {

    ERR_INIT

    const int precisionMin = 1;
    e = ErrorUtils::isAboveThreshold(
            precision,
            precisionMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_precision = precision;

    const double mzMaxMin = 100.0;
    e = ErrorUtils::isAboveThreshold(
            mzMax,
            mzMaxMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_mzMax = mzMaxMin;

    const int iterationsMaxMin = 1;
    e = ErrorUtils::isAboveThreshold(
            iterationsMax,
            iterationsMaxMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_iterationsMax = iterationsMax;

    const double stopTolMin = 0;
    e = ErrorUtils::isAboveThreshold(
            stopTolerance,
            stopTolMin,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;
    m_stopTolerance = stopTolerance;

    ERR_RETURN
}

namespace {

    Eigen::SparseMatrix<double> buildDeconvolveMatrix(
            const QList<QVector<MS2Ion>> &targetsMs2Ions,
            int precision,
            double mzMax
            ) {

        const int rows = MathUtils::hashDecimal(mzMax, precision);
        const int cols = targetsMs2Ions.size();

        Eigen::SparseMatrix<double> mat(rows, cols);

        int currentTargetIndex = 0;
        for (const QVector<MS2Ion> &targetIons : targetsMs2Ions) {

            for (const MS2Ion &ion : targetIons) {

                const int mzHashed = MathUtils::hashDecimal(ion.mz, precision);

                if (mzHashed >= rows || currentTargetIndex >= cols){
                    continue;
                }

                mat.insert(mzHashed, currentTargetIndex) = ion.intensity;

            }

            currentTargetIndex++;
        }

        return mat;
    }

}//namespace
Err TandemSpectraDeconvolvotron::deconvolveTandemSpectra(
        const ScanPoints &scanPoints,
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPredictions,
        QMap<PeptideStringWithMods, double> *pepSeqVsWeight
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isNotEmpty(tandemPredictions); ree;

    const Eigen::SparseMatrix<double> mat = buildDeconvolveMatrix(
            tandemPredictions.values(),
            m_precision,
            m_mzMax
            );

    Eigen::VectorX<double> vecScanPoints = EigenUtils::convertQPointFVecToEigen(
            scanPoints,
            m_precision,
            m_mzMax
            );

    vecScanPoints /= vecScanPoints.maxCoeff();

    Eigen::VectorX<double> x(mat.cols());
    x.setZero();

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> lscg;
    lscg.setMaxIterations(m_iterationsMax);

    lscg.compute(mat);
    x = lscg.solve(vecScanPoints);

    const QList<PeptideStringWithMods> &keys = tandemPredictions.keys();

    for (int i = 0; i < keys.size(); i++) {

        const PeptideStringWithMods &peptideSequenceChargeKey = keys.at(i);
        const double score = x.coeff(i);
        pepSeqVsWeight->insert(peptideSequenceChargeKey, score);
    }

    qDebug() << "#iterations:" << lscg.iterations() << "estimated error: " << lscg.error();

    ERR_RETURN
}
