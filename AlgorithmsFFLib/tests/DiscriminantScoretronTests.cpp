#include <QtTest>
#include <QCoreApplication>

#include "CandidateScores.h"
#include "DiscriminantScoretron.h"

class DiscriminantScoretronTest : public QObject {
Q_OBJECT

public:
    DiscriminantScoretronTest() = default;

    ~DiscriminantScoretronTest() override = default;

private slots:

    static void setDisciminantScoreForCandidatesTest();



};

void DiscriminantScoretronTest::setDisciminantScoreForCandidatesTest() {

    ERR_INIT

    QVector<QVector<float>> targetsData = {
            {1.0, 1.0, 1.0, 1.0, 5.0},
            {1.0, 1.0, 1.0, 1.0, 5.0},
            {1.0, 1.0, 1.0, 1.0, 5.0}
    };

    QVector<QVector<float>> decoyData = {
            {1.0, 1.0, 1.0, 1.0, 1.0},
            {1.0, 1.0, 1.0, 1.0, 1.0},
            {1.0, 1.0, 1.0, 1.0, 1.0}
    };

    const QStringList seqs = {
            "FLOPS",
            "CHAUNCY",
            "THOR"
    };

    QVector<TargetDecoyCandidatePair> targetDecoyPairsTargets = {
            TargetDecoyCandidatePair(),
            TargetDecoyCandidatePair(),
            TargetDecoyCandidatePair()
    };

    QVector<TargetDecoyCandidatePair> targetDecoyPairsDecoys = {
            TargetDecoyCandidatePair(),
            TargetDecoyCandidatePair(),
            TargetDecoyCandidatePair()
    };

    QVector<CandidateScores> candidateScoresVec;
    for (int i = 0; i < targetsData.size(); i++) {


        CandidateScores cs;
        cs.featuresArray = targetsData.at(i);
        cs.isDecoy = false;
        cs.targetDecoyCandidatePair = &targetDecoyPairsTargets[i];
        cs.targetDecoyCandidatePair->setCharge(2);
        cs.targetDecoyCandidatePair->setPeptideStringWithMods(seqs.at(i));
        cs.targetKey = "1234";

        CandidateScores csDecoy;
        csDecoy.featuresArray = decoyData.at(i);
        csDecoy.isDecoy = true;
        csDecoy.targetDecoyCandidatePair = &targetDecoyPairsDecoys[i];
        csDecoy.targetDecoyCandidatePair->setCharge(2);
        csDecoy.targetDecoyCandidatePair->setPeptideStringWithMods(seqs.at(i));
        csDecoy.targetKey = "1234";

        candidateScoresVec.push_back(cs);
        candidateScoresVec.push_back(csDecoy);
    }

    QVector<CandidateScores*> candidateScoresPtrsVec;
    std::transform(
            candidateScoresVec.begin(),
            candidateScoresVec.end(),
            std::back_inserter(candidateScoresPtrsVec),
            [](CandidateScores &c){return &c;}
            );

    e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
            false,
            false,
            &candidateScoresPtrsVec
            );
    QCOMPARE(e, eNoError);

    //TODO this is not really much of a test.  replace arrays w/ values that will yield a non-zero disc score.
    for (const CandidateScores *cs : candidateScoresPtrsVec) {
        QCOMPARE(cs->discriminantScore, 0);
    }

}


QTEST_MAIN(DiscriminantScoretronTest)

#include "DiscriminantScoretronTest.moc"
