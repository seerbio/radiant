#ifndef RADIANT_NEURAL_NET_FAMILY_FOLDS_H
#define RADIANT_NEURAL_NET_FAMILY_FOLDS_H

#include <QHash>
#include <QString>
#include <QVector>

#include <algorithm>

namespace NeuralNetFamilyFolds {

inline QString canonicalFamily(QString originSequence) {
    QString family;
    family.reserve(originSequence.size());
    bool inModification = false;
    for (const QChar character : originSequence) {
        if (character == QChar('(') || character == QChar('[')) {
            inModification = true;
        }
        else if (character == QChar(')') || character == QChar(']')) {
            inModification = false;
        }
        else if (!inModification && character != QChar('_')) {
            family.push_back(character == QChar('I') ? QChar('L') : character);
        }
    }
    return family;
}

// Assigns complete peptide families to deterministic, approximately balanced
// folds. The returned row assignments are stable for a given set of family
// names and keep every row from one family in the same fold.
inline bool buildFoldAssignments(
        const QVector<QString> &familyNames,
        int foldCount,
        QVector<int> *rowFolds
        ) {
    if (rowFolds == nullptr || familyNames.isEmpty() || foldCount < 2
        || foldCount > familyNames.size()) {
        return false;
    }

    QHash<QString, QVector<int>> rowsByFamily;
    for (int row = 0; row < familyNames.size(); ++row) {
        if (familyNames.at(row).isEmpty()) {
            return false;
        }
        rowsByFamily[familyNames.at(row)].push_back(row);
    }

    if (rowsByFamily.size() < foldCount) {
        return false;
    }

    struct FamilyRows {
        QString name;
        QVector<int> rows;
    };

    QVector<FamilyRows> families;
    families.reserve(rowsByFamily.size());
    for (auto it = rowsByFamily.cbegin(); it != rowsByFamily.cend(); ++it) {
        families.push_back({it.key(), it.value()});
    }

    std::sort(
        families.begin(),
        families.end(),
        [](const FamilyRows &left, const FamilyRows &right) {
            if (left.rows.size() != right.rows.size()) {
                return left.rows.size() > right.rows.size();
            }
            return left.name < right.name;
        }
    );

    QVector<int> foldSizes(foldCount, 0);
    QVector<int> familyFold(families.size(), -1);
    for (int family = 0; family < families.size(); ++family) {
        int selectedFold = 0;
        for (int fold = 1; fold < foldCount; ++fold) {
            if (foldSizes.at(fold) < foldSizes.at(selectedFold)) {
                selectedFold = fold;
            }
        }
        familyFold[family] = selectedFold;
        foldSizes[selectedFold] += families.at(family).rows.size();
    }

    for (int foldSize : foldSizes) {
        if (foldSize == 0) {
            return false;
        }
    }

    rowFolds->fill(-1, familyNames.size());
    for (int family = 0; family < families.size(); ++family) {
        for (int row : families.at(family).rows) {
            (*rowFolds)[row] = familyFold.at(family);
        }
    }

    return std::all_of(
        rowFolds->cbegin(),
        rowFolds->cend(),
        [](int fold) { return fold >= 0; }
    );
}

} // namespace NeuralNetFamilyFolds

#endif // RADIANT_NEURAL_NET_FAMILY_FOLDS_H
