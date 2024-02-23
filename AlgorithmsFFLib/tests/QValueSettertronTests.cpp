#include "QValueSettertron.h"

#include "CandidateScores.h"
#include "Error.h"

#include <QtTest>

class QValueSettertronTests : public QObject
{
    Q_OBJECT

public:
    QValueSettertronTests() = default;
    ~QValueSettertronTests() override = default;

private slots:

    static void setQValueForCandidatesTest();

};

void QValueSettertronTests::setQValueForCandidatesTest() {

    ERR_INIT

    TargetDecoyCandidatePair targetDecoyCandidatePair0;
    targetDecoyCandidatePair0.setPeptideStringWithMods(PeptideStringWithMods("PEPTIDEMURT"));
    targetDecoyCandidatePair0.setCharge(2);

    CandidateScores c0;
    c0.discriminantScore = 11;
    c0.classifierScore = 0.1;
    c0.isDecoy = false;
    c0.targetDecoyCandidatePair = &targetDecoyCandidatePair0;

    TargetDecoyCandidatePair targetDecoyCandidatePair1;
    targetDecoyCandidatePair1.setPeptideStringWithMods(PeptideStringWithMods("PEPTIDE"));
    targetDecoyCandidatePair1.setCharge(2);

    CandidateScores c1;
    c1.discriminantScore = 10;
    c1.classifierScore = 0.15;
    c1.isDecoy = false;
    c1.targetDecoyCandidatePair = &targetDecoyCandidatePair1;

    TargetDecoyCandidatePair targetDecoyCandidatePair2;
    targetDecoyCandidatePair2.setPeptideStringWithMods(PeptideStringWithMods("PEPTIDER"));
    targetDecoyCandidatePair2.setCharge(2);

    CandidateScores c2;
    c2.discriminantScore = 5;
    c2.classifierScore = 0.25;
    c2.isDecoy = false;
    c2.targetDecoyCandidatePair = &targetDecoyCandidatePair2;
    

    TargetDecoyCandidatePair targetDecoyCandidatePairDecoy;
    targetDecoyCandidatePair1.setPeptideStringWithMods(PeptideStringWithMods("PEPTIDECOY"));
    targetDecoyCandidatePair1.setCharge(2);

    CandidateScores cDecoy;
    cDecoy.discriminantScore = 1;
    cDecoy.classifierScore = 1.0;
    cDecoy.isDecoy = true;
    cDecoy.targetDecoyCandidatePair = &targetDecoyCandidatePairDecoy;

    TargetDecoyCandidatePair targetDecoyCandidatePair3;
    targetDecoyCandidatePair3.setPeptideStringWithMods(PeptideStringWithMods("PEPTIDEBUM"));
    targetDecoyCandidatePair3.setCharge(2);

    CandidateScores c3;
    c3.discriminantScore = 0.5;
    c3.classifierScore = 1.0;
    c3.isDecoy = false;
    c3.targetDecoyCandidatePair = &targetDecoyCandidatePair3;

    QVector<CandidateScores*> candidateScores = {
            &c0,
            &c1,
            &c2,
            &cDecoy,
            &c3
    };

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::DiscriminantScore,
            &candidateScores
            );
    QCOMPARE(e, eNoError);

    const QVector<float> expectedQValues = {0, 0, 0, 1, 0.25};
    QCOMPARE(candidateScores.size(), expectedQValues.size());

    for (int i = 0; i < expectedQValues.size(); i++) {
        QCOMPARE(candidateScores.at(i)->qValue, expectedQValues.at(i));
    }

}


QTEST_MAIN(QValueSettertronTests)

#include "QValueSettertronTests.moc"
