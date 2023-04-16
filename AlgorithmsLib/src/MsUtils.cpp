//
// Created by Drucifer on 4/27/2022.
//

#include "MsUtils.h"

#include "ErrorUtils.h"
#include "ParallelUtils.h"

namespace {
    const auto sortAscMz = [](const QPointF &l, const QPointF &r){return l.x() < r.x();};
}
ExtractPoints MsUtils::extractPointsFromPoints(
        const QVector<QPointF> &_points,
        const QVector<QPointF> &_pointsToExtract,
        double extractionPPM
        ) {

    ExtractPoints extractPointsOutput;
    extractPointsOutput.mzFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});
    extractPointsOutput.intensityFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});

    QVector<QPointF> pointsToExtract = _pointsToExtract;
    std::sort(pointsToExtract.begin(), pointsToExtract.end(), sortAscMz);

    QVector<QPointF> points = _points;
    std::sort(points.begin(), points.end(), sortAscMz);

    int currentExtractionIndex = 0;
    double extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
    double extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

    extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;

    double xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
    double xLo = extractionPointX - xPPM;
    double xHi = extractionPointX + xPPM;

    for (const QPointF &xPoint : points) {

        const double xVal = xPoint.x();

        while (xVal > xHi) {

            currentExtractionIndex++;

            if (currentExtractionIndex > _pointsToExtract.size() - 1) {
                break;
            }

            extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
            extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

            xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
            xLo = extractionPointX - xPPM;
            xHi = extractionPointX + xPPM;

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;
        }

        if (xLo <= xVal && xVal <= xHi) {

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex]
                    =  xPoint.y() > extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    ? QPointF(xVal , extractionPointX)
                    : extractPointsOutput.mzFoundVsSearched[currentExtractionIndex];

            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].rx() =  std::max(
                    xPoint.y(),
                    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    );
        }
    }

    return extractPointsOutput;
}

namespace {


     Err normalizeTheoPointsToActual(
            const ExtractPoints &extractPoints,
            double mzHiPassCutoff,
            QVector<QPointF> *normalizedPoints
            ) {


         ERR_INIT

         const QVector<QPointF> &mzVals = extractPoints.mzFoundVsSearched;
         const QVector<QPointF> &intzVals = extractPoints.intensityFoundVsSearched;

         const QPair<QVector<double>, QVector<double>> splitPointsMz = ParallelUtils::unZip(mzVals);
         const QPair<QVector<double>, QVector<double>> splitPointsIntz = ParallelUtils::unZip(intzVals);

         QVector<QPointF> founds;
         e = ParallelUtils::zip(
                 splitPointsMz.first,
                 splitPointsIntz.first,
                 &founds
         ); ree;

         QVector<QPointF> theos;
         e = ParallelUtils::zip(
                 splitPointsMz.second,
                 splitPointsIntz.second,
                 &theos
         ); ree;

         e = ErrorUtils::isEqual(
                 founds.size(),
                 theos.size()
         ); ree;

         QVector<double> sumFounds;
         QVector<double> sumTheos;
         for (int i = 0; i < founds.size(); ++i) {

             const QPointF &found = founds.at(i);
             const QPointF &theo = theos.at(i);

             if (found.x() > mzHiPassCutoff) {
                 sumFounds.push_back(found.y());
             }

             if (theo.x() > mzHiPassCutoff) {
                 sumTheos.push_back(theo.y());
             }
         }

         const double meanFounds = MathUtils::mean(sumFounds);
         const double meanTheos = MathUtils::mean(sumTheos);

         const double scaleFactor = MathUtils::tZero(meanTheos) ? 0.0 : meanFounds / meanTheos;

         QVector<QPointF> normalizedFoundPoints;
         for (int i = 0; i < splitPointsMz.first.size(); i++) {

             const double mz = splitPointsMz.first.at(i);
             const double intz = splitPointsIntz.second.at(i) * scaleFactor;
             normalizedFoundPoints.push_back({mz, intz});
         }

         *normalizedPoints = normalizedFoundPoints;

         ERR_RETURN
    }


}//namespace
Err MsUtils::buildDeletionPoints(
        const ExtractPoints &ep,
        double mzHiPassCutoff,
        QVector<QPointF> *delPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ep.mzFoundVsSearched); ree;
    e = ErrorUtils::isEqual(
            ep.intensityFoundVsSearched.size(),
            ep.mzFoundVsSearched.size()
            ); ree;

    QVector<QPointF> founds;
    e = normalizeTheoPointsToActual(
            ep,
            mzHiPassCutoff,
            &founds
            );

    ERR_RETURN
}


//Eigen::MatrixXd MsUtils::convertQPointsToMatrix(const QVector<QPointF> &points) {
//
//    Eigen::MatrixXd mat(points.size(), 2);
//    for (int i = 0; i < points.size(); i++) {
//        const QPointF &p = points.at(i);
//        mat.coeffRef(i, 0) = p.x();
//        mat.coeffRef(i, 1) = 1.0;
//    }
//
//    return mat;
//}
//
//void MsUtils::convertQPointsToSparseMatrix(
//        const QMap<int, QVector<QPointF>> &pointsMap,
//        double maxValue,
//        int precision,
//        int expectedPointsPerRow,
//        SparseMatrixD *output
//) {
//
//    const int buffer = 10;
//    const int maxValueHashed = MathUtils::hashDecimal(maxValue, precision);
//
//    const int rows = maxValueHashed + buffer;
//    const int cols = pointsMap.size();
//
//    SparseMatrixD mat(rows, cols);
//    mat.reserve(mat.rows() * expectedPointsPerRow);
//
//    for (auto it = pointsMap.begin(); it != pointsMap.end(); it++) {
//
//        const int index = it.key();
//        const QVector<QPointF> &points = it.value();
//
//        for (const QPointF &p : points) {
//
//            const int xValHashed = MathUtils::hashDecimal(p.x(), precision);
//
//            if (xValHashed > maxValueHashed) {
//                continue;
//            }
//
//            mat.insert(xValHashed, index) = p.y();
//        }
//
//    }
//
////    mat.conservativeResize(mat.rows(), mat.cols()); # TODO see if this helps w/ memory usuage.
//    *output = mat;
//}
//

//namespace {
//
//    struct ParallelSmoothLogicInput {
//        Eigen::MatrixXd kernelGauss;
//        Eigen::SparseVector<double> ogVec;
//        int index = -1;
//    };
//
//    struct ParallelSmoothLogicOutput {
//        Eigen::SparseVector<double> smoothedVec;
//        int index = -1;
//    };
//
//    ParallelSmoothLogicOutput parallelSmoothLogic(const ParallelSmoothLogicInput &input) {
//
//        const int kernelSize = static_cast<int>(input.kernelGauss.size());
//
//        const Eigen::SparseVector<double> newVec
//                = EigenKernelUtils::addPaddingToSparseVector(input.ogVec, kernelSize);
//
//        const Eigen::SparseVector<double> vv
//                = EigenKernelUtils::convolveVectorWithKernel(newVec, input.kernelGauss);
//
//        ParallelSmoothLogicOutput output;
//        output.index = input.index;
//        output.smoothedVec = vv;
//
//        return output;
//    }
//
//} //NAMESPACE
//Eigen::SparseMatrix<double> MsUtils::gaussianSmoothSparseMatrixColWise(
//        const Eigen::SparseMatrix<double> &mat,
//        int kernelLength,
//        double kernelStDev,
//        int expectedColObjects
//) {
//
//    const Eigen::MatrixXd kernelGauss
//            = EigenKernelUtils::buildGaussianFilter1D(kernelLength, kernelStDev);
//
//    QVector<ParallelSmoothLogicInput> matrixColumns;
//    for(int i = 0; i < mat.cols(); i++) {
//
//        ParallelSmoothLogicInput psli;
//        psli.kernelGauss = kernelGauss;
//        psli.ogVec = mat.col(i);
//        psli.index = i;
//
//        matrixColumns.push_back(psli);
//    }
//
//    QFuture<ParallelSmoothLogicOutput> future = QtConcurrent::mapped(
//            matrixColumns,
//            parallelSmoothLogic
//    );
//    future.waitForFinished();
//
//    Eigen::SparseMatrix<double> matSmoothed(mat.rows(), mat.cols());
//    matSmoothed.reserve(expectedColObjects * mat.cols());
//
//    for (const ParallelSmoothLogicOutput &op : future) {
//
//        const Eigen::SparseVector<double> &vec = op.smoothedVec;
//
//        for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it) {
//            const double &val = it.value();
//            const int &ogIndex = it.index();
//
//            matSmoothed.insert(ogIndex, op.index) = val;
//        }
//
//    }
//
//    return matSmoothed;
//}
//
//
//Eigen::SparseMatrix<double> MsUtils::gaussianSmoothSparseMatrixRowWise(
//        const Eigen::SparseMatrix<double> &_mat,
//        int kernelLength,
//        double kernelStDev,
//        int expectedRowObjects
//) {
//
//    Eigen::SparseMatrix<double, Eigen::RowMajor> mat = _mat;
//
//    const Eigen::MatrixXd kernelGauss
//            = EigenKernelUtils::buildGaussianFilter1D(kernelLength, kernelStDev);
//
//    QVector<ParallelSmoothLogicInput> matrixRows;
//    for(int i = 0; i < mat.rows(); i++) {
//
//        ParallelSmoothLogicInput psli;
//        psli.kernelGauss = kernelGauss;
//        psli.ogVec = mat.row(i);
//        psli.index = i;
//
//        matrixRows.push_back(psli);
//    }
//
//    QFuture<ParallelSmoothLogicOutput> future = QtConcurrent::mapped(
//            matrixRows,
//            parallelSmoothLogic
//    );
//    future.waitForFinished();
//
//    Eigen::SparseMatrix<double, Eigen::RowMajor> matSmoothed(mat.rows(), mat.cols());
//    matSmoothed.reserve(expectedRowObjects * mat.rows());
//
//    for (const ParallelSmoothLogicOutput &op : future) {
//
//        const Eigen::SparseVector<double> &vec = op.smoothedVec;
//
//        for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it) {
//            const double &val = it.value();
//            const int &ogIndex = it.index();
//            matSmoothed.insert(op.index, ogIndex) = val;
//        }
//
//    }
//
//    Eigen::SparseMatrix<double> matSmoothedColMajor = matSmoothed;
//    return matSmoothedColMajor;
//}
