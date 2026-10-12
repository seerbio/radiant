#include "NeuralNetFamilyFolds.h"

#include <QtTest>

class NeuralNetFamilyFoldsTests : public QObject
{
    Q_OBJECT

private slots:
    void canonicalizesEquivalentFamilies();
    void keepsFamiliesTogetherAndBalancesRows();
    void rejectsInvalidFoldRequests();
};

void NeuralNetFamilyFoldsTests::canonicalizesEquivalentFamilies() {
    QCOMPARE(
        NeuralNetFamilyFolds::canonicalFamily(QStringLiteral("AIC_D")),
        QStringLiteral("ALCD")
    );
    QCOMPARE(
        NeuralNetFamilyFolds::canonicalFamily(QStringLiteral("AL[+80]CD")),
        QStringLiteral("ALCD")
    );
    QVERIFY(
        NeuralNetFamilyFolds::canonicalFamily(QStringLiteral("PEPTIDE"))
        != NeuralNetFamilyFolds::canonicalFamily(QStringLiteral("PEPTIDK"))
    );
}

void NeuralNetFamilyFoldsTests::keepsFamiliesTogetherAndBalancesRows() {
    const QVector<QString> families{
        QStringLiteral("A"), QStringLiteral("A"), QStringLiteral("A"),
        QStringLiteral("B"),
        QStringLiteral("C"), QStringLiteral("C"),
        QStringLiteral("D"), QStringLiteral("D"),
        QStringLiteral("E"), QStringLiteral("F")
    };

    QVector<int> folds;
    QVERIFY(NeuralNetFamilyFolds::buildFoldAssignments(families, 3, &folds));
    QCOMPARE(folds.size(), families.size());

    QCOMPARE(folds.at(0), folds.at(1));
    QCOMPARE(folds.at(1), folds.at(2));
    QCOMPARE(folds.at(4), folds.at(5));
    QCOMPARE(folds.at(6), folds.at(7));

    QVector<int> foldCounts(3, 0);
    for (int fold : folds) {
        QVERIFY(fold >= 0 && fold < foldCounts.size());
        ++foldCounts[fold];
    }
    for (int count : foldCounts) {
        QVERIFY(count > 0);
    }

    QVector<int> repeatedFolds;
    QVERIFY(NeuralNetFamilyFolds::buildFoldAssignments(families, 3, &repeatedFolds));
    QCOMPARE(repeatedFolds, folds);
}

void NeuralNetFamilyFoldsTests::rejectsInvalidFoldRequests() {
    const QVector<QString> families{
        QStringLiteral("A"), QStringLiteral("A"), QStringLiteral("B")
    };
    QVector<int> folds;

    QVERIFY(!NeuralNetFamilyFolds::buildFoldAssignments(families, 1, &folds));
    QVERIFY(!NeuralNetFamilyFolds::buildFoldAssignments(families, 3, &folds));
    QVERIFY(!NeuralNetFamilyFolds::buildFoldAssignments(
        {QStringLiteral("A"), QString()}, 2, &folds));
    QVERIFY(!NeuralNetFamilyFolds::buildFoldAssignments(families, 2, nullptr));
}

QTEST_MAIN(NeuralNetFamilyFoldsTests)

#include "NeuralNetFamilyFoldsTests.moc"
