#include <QtTest>
#include <QCoreApplication>

#include "CandidateScores.h"
#include "Deconvolvotron.h"
#include "MsReaderParquet.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

#include "Eigen/Dense"
#include <unsupported/Eigen/NNLS>

#include <iostream>

class DeconvolvotronTests : public QObject {
Q_OBJECT

public:
    DeconvolvotronTests() = default;

    ~DeconvolvotronTests() override = default;

private slots:

    static void deconvoleTest();

};

void DeconvolvotronTests::deconvoleTest() {

    ERR_INIT

    const QVector<QPointF> scanPoints = {
        {100.01, 10.1},
        {200.01, 10},
        {300.01, 20},
        {400.01, 10},
        {500.01, 1}
    };

    const QVector<QPointF> peptide1 = {
        {100.01, 10},
        {200.01, 10},
        {300.01, 10}
    };

    const QVector<QPointF> peptide2 = {
        {300.01, 10},
        {400.01, 10},
    };

    const QVector<QPointF> peptide3 = {
        {310.01, 10},
        {410.01, 10},
    };

    CandidateScores cs1;
    CandidateScores cs2;
    CandidateScores cs3;

    QVector<QPair<CandidateScores*, QVector<QPointF>>> candScorePntrsVsPoints = {
        {&cs1, peptide1},
        {&cs2, peptide2},
        {&cs3, peptide3}
    };

    Deconvolvotron deconvolvotron;
    e = deconvolvotron.init(4, 10.0);
    QCOMPARE(e, eNoError);

    QVector<QPair<CandidateScores*, DeconvolvotronResult>> deconvolvotronResult;
    e = deconvolvotron.deconvolve(
        candScorePntrsVsPoints,
        scanPoints,
        &deconvolvotronResult
        );
    QCOMPARE(e, eNoError);

    QCOMPARE(deconvolvotronResult.size(), 2);

    //TODO fix pValue code

    qDebug() << deconvolvotronResult.at(0).second.discScore
             << deconvolvotronResult.at(0).second.pVal
             << deconvolvotronResult.at(0).second.pValFrameFtest;

    QCOMPARE(deconvolvotronResult.at(0).first, &cs1);
    QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.discScore, 0.5));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.pVal, 2.06302e-05));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.pValFrameFtest, 0.000492698));

    QCOMPARE(deconvolvotronResult.at(1).first, &cs2);
    QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.discScore, 0.5));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.pVal, 3.78698e-05));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.pValFrameFtest, 0.000492698));

}


QTEST_MAIN(DeconvolvotronTests)

#include "DeconvolvotronTests.moc"
