#include "QValueSettertron.h"

#include "CandidateScores.h"
#include "Error.h"

#include <QtTest>

#include <cmath>
#include <memory>
#include <vector>

namespace {

    class CandidateScorePairFixture {
    public:
        CandidateScorePairFixture(
            const QString &peptide,
            int charge,
            bool originalLibraryDecoy,
            double targetDiscriminantScore,
            double decoyDiscriminantScore,
            const MzTargetKey &targetKey
            )
        : targetDecoyCandidatePair(PeptideStringWithMods(peptide), 1.0f)
        {
            fragLibReaderRow.peptideSequenceChargeKey = peptide + "_" + QString::number(charge);
            fragLibReaderRow.isDecoy = originalLibraryDecoy ? 1 : 0;
            fragLibReaderRow.precursorCharge = charge;
            fragLibReaderRow.mass = 1000.0;
            fragLibReaderRow.iRT = 10.0;
            fragLibReaderRow.iM = 1.0;

            targetDecoyCandidatePair.setFragLibReaderRowPntr(&fragLibReaderRow);

            targetScores.targetDecoyCandidatePair = &targetDecoyCandidatePair;
            targetScores.isDecoy = false;
            targetScores.discriminantScore = targetDiscriminantScore;
            targetScores.classifierScore = 1.0 / targetDiscriminantScore;
            targetScores.targetKey = targetKey;

            decoyScores.targetDecoyCandidatePair = &targetDecoyCandidatePair;
            decoyScores.isDecoy = true;
            decoyScores.discriminantScore = decoyDiscriminantScore;
            decoyScores.classifierScore = 1.0 / decoyDiscriminantScore;
            decoyScores.targetKey = targetKey;
        }

        FragLibReaderRow fragLibReaderRow;
        TargetDecoyCandidatePair targetDecoyCandidatePair;
        CandidateScores targetScores;
        CandidateScores decoyScores;
    };

    std::unique_ptr<CandidateScorePairFixture> makeFixture(
        const QString &peptide,
        double targetDiscriminantScore,
        double decoyDiscriminantScore,
        bool originalLibraryDecoy = false,
        const MzTargetKey &targetKey = QStringLiteral("500000"),
        int charge = 2
        ) {

        return std::make_unique<CandidateScorePairFixture>(
            peptide,
            charge,
            originalLibraryDecoy,
            targetDiscriminantScore,
            decoyDiscriminantScore,
            targetKey
            );
    }

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> buildPairPointers(
        const std::vector<std::unique_ptr<CandidateScorePairFixture>> &fixtures
        ) {

        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs;
        pairs.reserve(static_cast<int>(fixtures.size()));
        for (const std::unique_ptr<CandidateScorePairFixture> &fixture : fixtures) {
            pairs.push_back({&fixture->targetScores, &fixture->decoyScores});
        }

        return pairs;
    }

    void compareDouble(double actual, double expected) {
        QVERIFY2(std::abs(actual - expected) < 1e-9, qPrintable(QStringLiteral("actual %1 expected %2").arg(actual).arg(expected)));
    }

}//namespace

class QValueSettertronTests : public QObject
{
    Q_OBJECT

public:
    QValueSettertronTests() = default;
    ~QValueSettertronTests() override = default;

private slots:

    static void setQValueForCandidatePairsUsesLegacyQValuesByDefaultTest();
    static void setQValueForCandidatePairsUsesMonotoneQValuesTest();
    static void setQValueForCandidatePairsCanStratifyMonotoneQValuesByTargetKeyTest();
    static void setQValueForCandidatePairsStratifiesMonotoneQValuesByChargeWithinTargetKeyTest();
    static void setQValueForCandidatePairsTreatsOriginalLibraryDecoysAsDecoysTest();

};

void QValueSettertronTests::setQValueForCandidatePairsUsesLegacyQValuesByDefaultTest() {

    ERR_INIT

    std::vector<std::unique_ptr<CandidateScorePairFixture>> fixtures;
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEA"), 100.0, 95.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEB"), 90.0, 10.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEC"), 80.0, 9.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDED"), 70.0, 8.0));

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs = buildPairPointers(fixtures);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs
        );
    QCOMPARE(e, eNoError);

    compareDouble(fixtures.at(0)->targetScores.qValue, 1.0);
    compareDouble(fixtures.at(1)->targetScores.qValue, 0.5);
    compareDouble(fixtures.at(2)->targetScores.qValue, 0.5);
    compareDouble(fixtures.at(3)->targetScores.qValue, 0.5);
    compareDouble(fixtures.at(3)->targetScores.decoyRatio, 0.5);

}

void QValueSettertronTests::setQValueForCandidatePairsUsesMonotoneQValuesTest() {

    ERR_INIT

    std::vector<std::unique_ptr<CandidateScorePairFixture>> fixtures;
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEA"), 100.0, 95.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEB"), 90.0, 10.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEC"), 80.0, 9.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDED"), 70.0, 8.0));

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs = buildPairPointers(fixtures);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true
        );
    QCOMPARE(e, eNoError);

    compareDouble(fixtures.at(0)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(1)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(2)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(3)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(3)->targetScores.decoyRatio, 0.0);

}

void QValueSettertronTests::setQValueForCandidatePairsCanStratifyMonotoneQValuesByTargetKeyTest() {

    ERR_INIT

    std::vector<std::unique_ptr<CandidateScorePairFixture>> fixtures;
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEA"), 100.0, 10.0, false, QStringLiteral("500000")));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEDECOY"), 80.0, 7.0, true, QStringLiteral("500000")));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEB"), 70.0, 1.0, false, QStringLiteral("600000")));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEC"), 60.0, 1.0, false, QStringLiteral("600000")));

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs = buildPairPointers(fixtures);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true,
        false
        );
    QCOMPARE(e, eNoError);
    compareDouble(fixtures.at(2)->targetScores.qValue, 1.0 / 3.0);
    compareDouble(fixtures.at(3)->targetScores.qValue, 1.0 / 3.0);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true,
        true
        );
    QCOMPARE(e, eNoError);

    compareDouble(fixtures.at(0)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(1)->targetScores.qValue, 1.0);
    compareDouble(fixtures.at(2)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(3)->targetScores.qValue, 0.0);
}

void QValueSettertronTests::setQValueForCandidatePairsStratifiesMonotoneQValuesByChargeWithinTargetKeyTest() {

    ERR_INIT

    const MzTargetKey targetKey = QStringLiteral("500000");
    std::vector<std::unique_ptr<CandidateScorePairFixture>> fixtures;
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEA"), 100.0, 10.0, false, targetKey, 2));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEDECOY"), 90.0, 7.0, true, targetKey, 2));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEB"), 80.0, 1.0, false, targetKey, 3));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEC"), 70.0, 1.0, false, targetKey, 3));

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs = buildPairPointers(fixtures);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true,
        false
        );
    QCOMPARE(e, eNoError);
    QVERIFY(fixtures.at(2)->targetScores.qValue > 0.0);
    QVERIFY(fixtures.at(3)->targetScores.qValue > 0.0);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true,
        true
        );
    QCOMPARE(e, eNoError);

    compareDouble(fixtures.at(0)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(1)->targetScores.qValue, 1.0);
    compareDouble(fixtures.at(2)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(3)->targetScores.qValue, 0.0);
}

void QValueSettertronTests::setQValueForCandidatePairsTreatsOriginalLibraryDecoysAsDecoysTest() {

    ERR_INIT

    std::vector<std::unique_ptr<CandidateScorePairFixture>> fixtures;
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEA"), 100.0, 10.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEB"), 90.0, 9.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEC"), 80.0, 8.0));
    fixtures.push_back(makeFixture(QStringLiteral("PEPTIDEDECOY"), 95.0, 7.0, true));

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairs = buildPairPointers(fixtures);

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        &pairs,
        true
        );
    QCOMPARE(e, eNoError);

    compareDouble(fixtures.at(0)->targetScores.qValue, 0.0);
    compareDouble(fixtures.at(1)->targetScores.qValue, 1.0 / 3.0);
    compareDouble(fixtures.at(2)->targetScores.qValue, 1.0 / 3.0);
    compareDouble(fixtures.at(3)->targetScores.qValue, 1.0);
    compareDouble(fixtures.at(2)->targetScores.decoyRatio, 1.0 / 3.0);

}

QTEST_MAIN(QValueSettertronTests)

#include "QValueSettertronTests.moc"
