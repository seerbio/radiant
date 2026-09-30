#include "IdLevelQValueAnnotator.h"

#include "CandidateScores.h"
#include "FragLibReaderRow.h"
#include "TargetDecoyCandidatePair.h"

#include <QtTest/QtTest>

class IdLevelQValueAnnotatorTests : public QObject
{
    Q_OBJECT

private Q_SLOTS:
    static void annotatesCollapsedLevelsAndBestRows();
};

namespace {

    struct CandidateFixture {
        FragLibReaderRow fragLibReaderRow;
        TargetDecoyCandidatePair targetDecoyCandidatePair;
        CandidateScores candidateScores;
    };

    CandidateFixture makeCandidate(
        const QString &peptide,
        int charge,
        const QString &proteinGroup,
        bool isSyntheticDecoy,
        bool isLibraryDecoy,
        double classifierScore
        ) {

        CandidateFixture fixture;
        fixture.fragLibReaderRow.precursorCharge = charge;
        fixture.fragLibReaderRow.proteinGroups = proteinGroup;
        fixture.fragLibReaderRow.isDecoy = isLibraryDecoy;
        fixture.targetDecoyCandidatePair = TargetDecoyCandidatePair(PeptideStringWithMods(peptide), 0.0f);
        fixture.targetDecoyCandidatePair.setFragLibReaderRowPntr(&fixture.fragLibReaderRow);
        fixture.candidateScores.targetDecoyCandidatePair = &fixture.targetDecoyCandidatePair;
        fixture.candidateScores.proteinGroup = proteinGroup;
        fixture.candidateScores.isDecoy = isSyntheticDecoy;
        fixture.candidateScores.classifierScore = classifierScore;
        fixture.candidateScores.discriminantScore = 1.0 / classifierScore;
        return fixture;
    }

} // namespace

void IdLevelQValueAnnotatorTests::annotatesCollapsedLevelsAndBestRows() {

    CandidateFixture targetBest = makeCandidate("PEPTIDEK", 2, "P00001", false, false, 0.01);
    CandidateFixture targetDuplicate = makeCandidate("PEPTIDEK", 2, "P00001", false, false, 0.02);
    CandidateFixture targetWorseThanDecoy = makeCandidate("OTHERPEP", 2, "P00002", false, false, 0.90);
    CandidateFixture syntheticDecoy = makeCandidate("DECOYPEP", 2, "P00002", true, false, 0.80);

    QVector<CandidateScores*> candidates = {
        &targetBest.candidateScores,
        &targetDuplicate.candidateScores,
        &targetWorseThanDecoy.candidateScores,
        &syntheticDecoy.candidateScores
    };

    Err e = IdLevelQValueAnnotator::annotate(&candidates);
    QCOMPARE(e, eNoError);

    QCOMPARE(targetBest.candidateScores.isBestPrecursorCandidate, 1);
    QCOMPARE(targetBest.candidateScores.isBestPeptideCandidate, 1);
    QCOMPARE(targetBest.candidateScores.isBestProteinCandidate, 1);

    QCOMPARE(targetDuplicate.candidateScores.isBestPrecursorCandidate, 0);
    QCOMPARE(targetDuplicate.candidateScores.isBestPeptideCandidate, 0);
    QCOMPARE(targetDuplicate.candidateScores.isBestProteinCandidate, 0);

    QVERIFY(targetBest.candidateScores.precursorQValue <= targetDuplicate.candidateScores.precursorQValue);
    QCOMPARE(targetBest.candidateScores.proteinQValue, 0.0);
    QCOMPARE(targetWorseThanDecoy.candidateScores.proteinQValue, 0.5);
}

QTEST_MAIN(IdLevelQValueAnnotatorTests)

#include "IdLevelQValueAnnotatorTests.moc"
