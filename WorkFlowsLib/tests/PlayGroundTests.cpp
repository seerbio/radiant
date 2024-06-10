//
// Created by anichols on 11/07/2021.
//

#include "../../UtilsLib/src/ErrorUtils.h"
#include "../../FileReadersLib/src/MsReaderBase.h"

#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>


#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "CandidateScores.h"
#include "Deconvolvotron.h"
#include "GlobalSettings.h"
#include "FragLibReader.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePairManager.h"


class PlayGroundTests : public QObject
{
    Q_OBJECT

public:
    PlayGroundTests() = default;
    ~PlayGroundTests() override = default;

private Q_SLOTS:

void testme();

};

void PlayGroundTests::testme() {

    ERR_INIT

    const ScanNumber scanNumber = 11114;
    const QString targetKey = "625034";

    const QString fragLibFilepath = QStringLiteral("/home/anichols/Desktop/Data/Libraries/predicted_lib.speclib.fragLibFF");
    const QString msDataFilepath = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML.prqFF");

    MsReaderPointerAcc msReader;
    e = msReader.openFile(msDataFilepath);
    QCOMPARE(e, eNoError);

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(fragLibFilepath, 400, 4000, &fragLibReaderRows);
    QCOMPARE(e, eNoError);

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(PythiaParameterReader::genericPythiaParametersForTests(), &fragLibReaderRows);
    QCOMPARE(e, eNoError);

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs;
    e = targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&targetDecoyCandidatePairsPntrs);
    QCOMPARE(e, eNoError);

    const QStringList peps = {
        "SPYTVTVGQACNPSACR",
        "FLPSELRDEH"
    };

    TargetDecoyCandidatePair* pep1;
    TargetDecoyCandidatePair* pep2;
    for (TargetDecoyCandidatePair* pr : targetDecoyCandidatePairsPntrs) {
        if (pr->peptideStringWithMods() == peps.front() && pr->charge() == 3) {
            pep1 = pr;
        }
        if (pr->peptideStringWithMods() == peps.back() && pr->charge() == 2) {
            pep2 = pr;
        }
    }

    CandidateScores cs1;
    cs1.targetDecoyCandidatePair = pep1;
    cs1.isDecoy = false;
    CandidateScores cs2;
    cs2.targetDecoyCandidatePair = pep2;
    cs2.isDecoy = false;
    CandidateScores cs3;
    cs3.targetDecoyCandidatePair = pep2;
    cs3.isDecoy = true;
    // cs3.targetDecoyCandidatePair->ms2IonsTarget() = cs3.targetDecoyCandidatePair->ms2IonsDecoy();

    QPair<Err, ScanPoints*> scanPoints = msReader.ptr->getScanPoints(scanNumber);
    QCOMPARE(e, eNoError);

    QVector<QPointF> bVecPoints;
    std::transform(
        scanPoints.second->begin(),
        scanPoints.second->end(),
        std::back_inserter(bVecPoints),
        [](const ScanPoint &sp){return QPointF(sp.x(), sp.y());}
        );

    Ms2IonFraggertronManager ms2IonFraggertronManager;
    e = ms2IonFraggertronManager.initTesting({&cs1, &cs2});
    QCOMPARE(e, eNoError);

    QHash<CandidateScores*, QVector<MS2Ion>> candScoresVsMs2Ions;

    const float ppmTol = 10.0;
    for (const ScanPoint &sp : *scanPoints.second) {

        const float massTol = MathUtils::calculatePPM(sp.x(), ppmTol);

        QVector<QPair<MS2Ion, CandidateScores*>> ms2IonsVsCandidateScoresPntrses;
        e = ms2IonFraggertronManager.extractMs2Points(
            sp.x() - massTol,
            sp.x() + massTol,
            45.0,
            48.0,
            &ms2IonsVsCandidateScoresPntrses
            );

        for (const QPair<MS2Ion, CandidateScores*> &pr : ms2IonsVsCandidateScoresPntrses) {
            candScoresVsMs2Ions[pr.second].push_back(pr.first);
        }
    }

    QVector<QPair<CandidateScores*, QVector<QPointF>>> aMatrixPoints;
    for (auto it = candScoresVsMs2Ions.begin(); it != candScoresVsMs2Ions.end(); ++it) {
        QVector<QPointF> pnts;
        std::transform(
            it.value().begin(),
            it.value().end(),
            std::back_inserter(pnts),
            [](const MS2Ion &ms2Ion){return QPointF(ms2Ion.mz, ms2Ion.intensity);}
            );
        aMatrixPoints.push_back({it.key(), pnts});
    }

    Deconvolvotron deconvolvotron;
    e = deconvolvotron.init(1, 20.0);
    QCOMPARE(e, eNoError);

    QVector<QPair<CandidateScores*, DeconvolvotronResult>> candScoresVsScore;
    e = deconvolvotron.deconvolve(aMatrixPoints, bVecPoints, &candScoresVsScore);
    for (const QPair<CandidateScores*, DeconvolvotronResult> &res : candScoresVsScore) {
        qDebug() << res.second.discScore << res.second.pVal << res.second.pValFrameFtest;
    }

    // QMap<int, QVector<double>> hashedXValsVsMzValsGrouped;
    // e = deconvolvotron.buildXValsSet(aMatrixPoints, &hashedXValsVsMzValsGrouped);
    // QCOMPARE(e, eNoError);
    //
    // QVector<QVector<double>> aMat;
    // QVector<double> bVec;
    // e = Deconvolvotron::buildAMatrixAndBVecTestAccess(
    //     aMatrixPoints,
    //     bVecPoints,
    //     hashedXValsVsMzValsGrouped,
    //     20,
    //     &aMat,
    //     &bVec
    //     );
    // QCOMPARE(e, eNoError);

}

QTEST_MAIN(PlayGroundTests)
#include "PlayGroundTests.moc"
