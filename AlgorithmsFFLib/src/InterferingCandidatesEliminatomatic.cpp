//
// Created by anichols on 1/9/24.
//

#include "InterferingCandidatesEliminatomatic.h"

#include "CandidateScores.h"
#include "EigenUtils.h"
#include "GlobalSettings.h"

Err InterferingCandidatesEliminatomatic::removeInterferingCandidates(
        int ionsSharedToReject,
        double mzMinMS2,
        double mzMaxMS2,
        QVector<CandidateScores*> *candidates
        ) {

    ERR_INIT

    QMap<ScanNumber, QVector<CandidateScores*>> scanNumberVsCandidateScoresPtr;
    for (int i = 0; i < candidates->size(); i++) {
        CandidateScores *cs = candidates->at(i);

        if (cs == nullptr) {
            continue;
        }

        scanNumberVsCandidateScoresPtr[cs->scanNumber].push_back(cs);
    }

    QVector<CandidateScores*> filteredCandidates;
    for (auto it = scanNumberVsCandidateScoresPtr.begin(); it != scanNumberVsCandidateScoresPtr.end(); it++) {

        const ScanNumber scanNumber = it.key();
        QVector<CandidateScores*> &scanNumberCandidates = it.value();
        const int foundCandidatesCountInScan = scanNumberCandidates.size();

        std::sort(
                scanNumberCandidates.rbegin(),
                scanNumberCandidates.rend(),
                [](CandidateScores *l, CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
        );

        if (foundCandidatesCountInScan == 1) {
            filteredCandidates.push_back(scanNumberCandidates.front());
            continue;
        }

        const int cols = static_cast<int>(std::round(mzMaxMS2 - mzMinMS2));
        const int rows = foundCandidatesCountInScan;
        Eigen::MatrixX<float> mat(rows, cols);
        mat.setZero();


        for (int row = 0; row < rows; row++) {

            const CandidateScores* cs = scanNumberCandidates.at(row);

            constexpr int top6MzValsFound = 12;
            const QVector<float> top6MzValsFoundArr = cs->featuresArray.mid(MzFoundMean1, top6MzValsFound);
            for (float mz : top6MzValsFoundArr) {

                const int col = static_cast<int>(std::round(mz));

                if (col < 0 || col >= cols) {
                    continue;
                }
                mat.coeffRef(row, col) = 1.0f;
            }
        }

        QSet<int> excludeIndexes;
        for (int row = 0; row < rows - 1; row++) {

            const Eigen::VectorX<float> &vecBase = mat.row(row);

            int rowNext = row + 1;
            const Eigen::MatrixX<float> &matNext = mat.block(rowNext, 0, rows - rowNext, cols);

            const Eigen::VectorX<float> dotProds =  matNext * vecBase;

            if (static_cast<int>(dotProds.maxCoeff()) < ionsSharedToReject) {
                continue;
            }

            for (int i = 0; i < dotProds.size(); i++) {
                const float sharedIons = dotProds.coeff(i);
                if (static_cast<int>(sharedIons) >= ionsSharedToReject) {
                    excludeIndexes.insert(row + i + 1);
                }
            }
        }

        for (int i = 0; i < scanNumberCandidates.size(); i++) {
            CandidateScores *cs = scanNumberCandidates.at(i);
            if (excludeIndexes.contains(i)) {
                continue;
            }
            filteredCandidates.push_back(cs);
        }
    }

    candidates->clear();
    *candidates = filteredCandidates;

    ERR_RETURN

}
