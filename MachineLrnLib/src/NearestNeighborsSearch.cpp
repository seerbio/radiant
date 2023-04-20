//
// Created by anichols on 4/18/23.
//

#include "NearestNeighborsSearch.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "EigenMatrixDiskUtils.h"
#include "GlobalSettings.h"

#include "../ThirdPartyLibs/nanoflann.hpp" //TODO find a better way to do this.

#include <Eigen/Dense>

#include <QElapsedTimer>
#include <QFile>
#include <QDataStream>

class Q_DECL_HIDDEN NearestNeighborsSearch::Private
{
public:

    using KDTree = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd>;

    Private();
    ~Private();

    Err init(
            const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
            int maxTreeLeafSize
            );

    Err kNearestNeighborsSearch(
            const QVector<QVector<double>> &searchPointCoors,
            int k,
            QVector<NNSearchResult> *result
    );

    Err radiusSearch(
            const QVector<QVector<double>> &searchPointCoors,
            double searchRadiusSquared,
            QVector<NNSearchResult> *result
    );

    [[nodiscard]] int treeSize() const;

    Err writeNearestNeighbors(
            const QString &msDataFilePath,
            QString *calibrationMatFilePath,
            QString *calibarationCalFilePath
            );

    Err readNearestNeighbors(
            const QString &calFilePath,
            const QString &matFilePath
            );

private:

    KDTree *m_kdTree;
    Eigen::MatrixX<double> *m_mat;
    int m_maxTreeLeafSize;
    int m_columnCount;

    QVector<double> m_pairValues;

};

NearestNeighborsSearch::Private::Private()
: m_maxTreeLeafSize(30)
, m_columnCount(0)
{}

NearestNeighborsSearch::Private::~Private() {
//    delete m_kdTree; //TODO, make sure this is not a memory leak
//    delete m_mat; //TODO, make sure this is not a memory leak
//    It appears d_ptr which is a qscopedpointer is deleting this so this is not
//    necessary.
}

namespace {

     Err extractColumnFromTreePoints(
            const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
            int columnIdx,
            int treeColumnCount,
            Eigen::VectorX<double> *eVec
                    ) {

         ERR_INIT


         const int rowCount = valuesVsTreePoints.size();

         e = ErrorUtils::isBelowThreshold(
                 columnIdx,
                 treeColumnCount,
                 ErrorUtilsParam::IncludeThreshold
                 ); ree;

         QVector<double> vec(rowCount);
         vec.reserve(rowCount);
         for (int i = 0; i < rowCount; i++) {

             e = ErrorUtils::isEqual(
                     valuesVsTreePoints[i].second.size(),
                     treeColumnCount
                     ); ree;

             vec[i] = valuesVsTreePoints[i].second[columnIdx];
         }

         *eVec = EigenUtils::convertQVectorToEigenVector(vec);
         ERR_RETURN
    }

    Err convertTreePointsToEigenMatrix(
            const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
            int treeColumnCount,
            Eigen::MatrixX<double> *mat
            ) {

         ERR_INIT

        e = ErrorUtils::isNotEmpty(valuesVsTreePoints); ree;

         const int columnCount = valuesVsTreePoints.front().second.size();
         const int rowCount = valuesVsTreePoints.size();

         mat->resize(rowCount, columnCount);

         for (int colIdx = 0; colIdx < columnCount; colIdx++) {

             Eigen::VectorX<double> eVec;
             e = extractColumnFromTreePoints(
                     valuesVsTreePoints,
                     colIdx,
                     treeColumnCount,
                     &eVec
                     ); ree;

             mat->col(colIdx) = eVec;
         }

         ERR_RETURN
    }

    Err extractValuesFromPairs(
            const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
            QVector<double> *vals
            ) {

         ERR_INIT

         e = ErrorUtils::isNotEmpty(valuesVsTreePoints); ree;

         vals->reserve(valuesVsTreePoints.size());
         for(const QPair<double, QVector<double>> &pr : valuesVsTreePoints) {
             vals->push_back(pr.first);
         }

         ERR_RETURN
     }

}//namespace
Err NearestNeighborsSearch::Private::init(
        const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
        int maxTreeLeafSize
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(valuesVsTreePoints); ree;
    m_maxTreeLeafSize = maxTreeLeafSize;

    m_columnCount = valuesVsTreePoints.front().second.size();
    e = ErrorUtils::isTrue(m_columnCount > 0); ree;

    m_mat = new Eigen::MatrixX<double>();

    e = convertTreePointsToEigenMatrix(
            valuesVsTreePoints,
            m_columnCount,
            m_mat
            ); ree;

    e = extractValuesFromPairs(
            valuesVsTreePoints,
            &m_pairValues
            ); ree;

    m_kdTree = new KDTree(
            m_columnCount,
            *m_mat,
            m_maxTreeLeafSize
            );

    const int treePointCount = static_cast<int>(m_kdTree->kdtree_get_point_count());
    e = ErrorUtils::isNotEqual(
            treePointCount,
            0
            ); ree;

    qDebug() << treePointCount << "Loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}

Err NearestNeighborsSearch::Private::kNearestNeighborsSearch(
        const QVector<QVector<double>> &searchPointCoors,
        int k,
        QVector<NNSearchResult> *result
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(searchPointCoors); ree;
    e = ErrorUtils::isNotEqual(
            static_cast<int>(m_kdTree->kdtree_get_point_count()),
            0
            ); ree; //TODO, figure out a better way to check this as this value is not
                        // initialized w/ nanoflann.

    result->reserve(searchPointCoors.size());

    const auto inserterLogic = [&](const long idx){
        return m_pairValues.at(static_cast<int>(idx));
    };

    for (const QVector<double> &coor : searchPointCoors) {

        std::vector<double> queryPt = coor.toStdVector();

        std::vector<long> retIndex(k);
        std::vector<double> outDistSqr(k);

        std::vector<std::pair<Eigen::Index, double>> matches;

        const size_t resultsSize = m_kdTree->index->knnSearch(
            queryPt.data(),
            k,
            retIndex.data(),
            outDistSqr.data()
            );

        retIndex.resize(resultsSize);
        outDistSqr.resize(resultsSize);

        std::vector<double> vals;
        vals.reserve(resultsSize);

        std::transform(
                retIndex.begin(),
                retIndex.end(),
                std::back_inserter(vals),
                inserterLogic
                );

        NNSearchResult nnSearchResult(
                coor,
                retIndex,
                outDistSqr,
                vals
                );

        result->push_back(nnSearchResult);
    }

    ERR_RETURN
}

Err NearestNeighborsSearch::Private::radiusSearch(
        const QVector<QVector<double>> &searchPointCoors,
        double searchRadiusSquared,
        QVector<NNSearchResult> *result
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(searchPointCoors); ree;
    e = ErrorUtils::isNotEqual(
            static_cast<int>(m_kdTree->kdtree_get_point_count()),
            0
    ); ree;

    nanoflann::SearchParams params;

    result->reserve(searchPointCoors.size());

    for (const QVector<double> &coor : searchPointCoors) {

        std::vector<std::pair<Eigen::Index, double>> foundPoints;

        const size_t resultsSize
                = m_kdTree->index->radiusSearch(
                        coor.data(),
                        searchRadiusSquared,
                        foundPoints,
                        params
                );

        NNSearchResult nnSearchResult;
        nnSearchResult.values.reserve(foundPoints.size());
        nnSearchResult.indexes.reserve(foundPoints.size());
        nnSearchResult.distancesSquared.reserve(foundPoints.size());

        for (const std::pair<Eigen::Index, double> &pr : foundPoints) {
            nnSearchResult.indexes.push_back(pr.first);
            nnSearchResult.values.push_back(m_pairValues.at(static_cast<int>(nnSearchResult.indexes.back())));
            nnSearchResult.distancesSquared.push_back(pr.second);
        }

        result->push_back(nnSearchResult);
    }

    ERR_RETURN
}

int NearestNeighborsSearch::Private::treeSize() const {
    return static_cast<int>(m_kdTree->kdtree_get_point_count());
}

Err NearestNeighborsSearch::Private::writeNearestNeighbors(
        const QString &msDataFilePath,
        QString *calibrationMatFilePath,
        QString *calibarationCalFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_pairValues); ree;
    e = ErrorUtils::isNotEqual(static_cast<int>(m_mat->rows()), 0); ree

    const QString matrixFilePath = msDataFilePath + S_GLOBAL_SETTINGS.DOT_MAT;
    e = EigenMatrixDiskUtils::saveEigenMatrix(
            *m_mat,
            matrixFilePath
            ); ree;

    qDebug() << "Mat file saved to" << matrixFilePath;


    const QString calFilePath = msDataFilePath + S_GLOBAL_SETTINGS.DOT_CAL;
    QFile writeFile(calFilePath);
    if (!writeFile.open(QIODevice::WriteOnly)) {
        qDebug() << "Could not write to file:" << calFilePath
                 << "Error string:" << writeFile.errorString();
        return Error::eFileError;
    }

    QDataStream out(&writeFile);
    out.setVersion(QDataStream::Qt_5_12);

    out << m_pairValues;
    out << m_maxTreeLeafSize;

    writeFile.close();

    qDebug() << "Cal file saved to" << calFilePath;

    *calibrationMatFilePath = matrixFilePath;
    *calibarationCalFilePath = calFilePath;

    ERR_RETURN
}

Err NearestNeighborsSearch::Private::readNearestNeighbors(
        const QString &calFilePath,
        const QString &matFilePath
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    m_mat = new Eigen::MatrixX<double>();

    e = EigenMatrixDiskUtils::loadEigenMatrix(
            matFilePath,
            m_mat
            ); ree;

    e = ErrorUtils::isFalse(static_cast<int>(m_mat->rows(), 0)); ree;

    QFile readFile(calFilePath);
    QDataStream in(&readFile);

    in.setVersion(QDataStream::Qt_5_12);

    if (!readFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Could not read the file:" << calFilePath
                 << "Error string:" << readFile.errorString();
        return Error::eFileError;
    }

    in >> m_pairValues;
    in >> m_maxTreeLeafSize;

    e = ErrorUtils::isNotEmpty(m_pairValues); ree;
    e = ErrorUtils::isNotEqual(m_maxTreeLeafSize, 0); ree;

    m_columnCount = static_cast<int>(m_mat->cols());

    m_kdTree = new KDTree(
            m_columnCount,
            *m_mat,
            m_maxTreeLeafSize
    );

    const int treePointCount = static_cast<int>(m_kdTree->kdtree_get_point_count());
    e = ErrorUtils::isNotEqual(
            treePointCount,
            0
    ); ree;

    qDebug() << treePointCount << "Loaded in" << et.elapsed() << "mSec";

    ERR_RETURN
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

NearestNeighborsSearch::NearestNeighborsSearch() : d_ptr(new Private()) {}

NearestNeighborsSearch::~NearestNeighborsSearch() {}

Err NearestNeighborsSearch::init(
        const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
        int maxTreeLeafSize
        ) {
    ERR_INIT

    e = d_ptr->init(
            valuesVsTreePoints,
            maxTreeLeafSize
            ); ree;

    ERR_RETURN
}

Err NearestNeighborsSearch::kNearestNeighborsSearch(
        const QVector<QVector<double>> &searchPointCoors,
        int k,
        QVector<NNSearchResult> *result
        ) {

    ERR_INIT

    e = d_ptr->kNearestNeighborsSearch(
            searchPointCoors,
            k,
            result
            ); ree;

    ERR_RETURN
}

Err NearestNeighborsSearch::radiusSearch(
        const QVector<QVector<double>> &searchPointCoors,
        double searchRadiusSquared,
        QVector<NNSearchResult> *result
        ) {
    ERR_INIT

    e = d_ptr->radiusSearch(
            searchPointCoors,
            searchRadiusSquared,
            result
            ); ree;

    ERR_RETURN
}

int NearestNeighborsSearch::kdTreeSize() const {
    return d_ptr->treeSize();
}

Err NearestNeighborsSearch::writeNearestNeighbors(
        const QString &msDataFilePath,
        QString *calibrationMatFilePath,
        QString *calibarationCalFilePath
        ) {

    ERR_INIT

    e = d_ptr->writeNearestNeighbors(
            msDataFilePath,
            calibrationMatFilePath,
            calibarationCalFilePath
            ); ree;

    ERR_RETURN
}

Err NearestNeighborsSearch::readNearestNeighbors(
        const QString &calFilePath,
        const QString &matFilePath
        ) {

    ERR_INIT
    e = d_ptr->readNearestNeighbors(
            calFilePath,
            matFilePath
            ); ree;
    ERR_RETURN
}

