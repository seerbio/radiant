//
// Created by anichols on 11/07/2021.
//

#include "EigenSparseUtils.h"


void  EigenSparseUtils::sortApexesHiLoValue(QVector<SparseMatrixPoint> *apexes)
{
    const auto sortLogic = [](const SparseMatrixPoint &l, const SparseMatrixPoint &r){
        return l.value < r.value;
    };
    std::sort(apexes->rbegin(), apexes->rend(), sortLogic);
}

Eigen::SparseMatrix<double> EigenSparseUtils::loadFrameToSparseMatrixColMajor(
        const QMap<int, QVector<QPointF>> &frame,
        int precision,
        double maxRowValue
        ) {

    const int maxRowValueHashed = MathUtils::hashDecimal(maxRowValue, precision);

    Eigen::SparseMatrix<double> mat(maxRowValueHashed, frame.size());

    for (auto it = frame.begin(); it != frame.end(); it++) {

        const int col =it.key();
        const QVector<QPointF> &points = it.value();

        for (const QPointF &p : points) {

            const double x = p.x();
            const int xHashed = MathUtils::hashDecimal(x, precision);
            const double y = p.y();

            if (xHashed < 0 || xHashed >= maxRowValueHashed || col >= frame.size()) {
                continue;
            }

            mat.coeffRef(xHashed, col) += y;
        }

    }

    return mat;
}

Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseUtils::loadFrameToSparseMatrixRowMajor(
        const QMap<int, QVector<QPointF>> &frame,
        int precision,
        double maxRowValue
        ) {

    const Eigen::SparseMatrix<double, Eigen::ColMajor> matColMaj = loadFrameToSparseMatrixColMajor(
            frame,
            precision,
            maxRowValue
            );

    return {matColMaj};
}



namespace {

    void sortFramePointsAsc(QMap<int, QVector<QPointF>> *points) {

        for (auto it = points->begin(); it != points->end(); it++) {

            QVector<QPointF> &p = it.value();

            std::sort(p.begin(), p.end(), [](const QPointF &l, const QPointF &r){
                return l.x() < l.x();
            });
        }

    }

}
QMap<int, QVector<QPointF>> EigenSparseUtils::loadSparseMatrixToFrame(
        const Eigen::SparseMatrix<double> &mat,
        int precision
        ) {

    QMap<int, QVector<QPointF>> frame;

    for (int i = 0; i < mat.outerSize(); ++i) {

        for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it){

            const int key = static_cast<int>(it.col());
            const auto x = MathUtils::unHashDecimal<double>(
                    static_cast<int>(it.row()),
                    precision
                    );
            const double val = it.value();
            frame[key].push_back({x, val});
        }
    }

//    sortFramePointsAsc(&frame);

    return frame;
}
