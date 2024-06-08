#include <QtTest>
#include <QCoreApplication>

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
        {100.01, 10},
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

    QVector<QPair<PeptideString, QVector<QPointF>>> pepStrVsPoints = {
        {"Pep1", peptide1},
        {"Pep2", peptide2},
        // {"Pep3", peptide3}
    };

    Deconvolvotron deconvolvotron;
    e = deconvolvotron.init(4);
    QCOMPARE(e, eNoError);

    QVector<QPair<IdStr, DeconvolvotronResult>> deconvolvotronResult;
    e = deconvolvotron.deconvolve(
        pepStrVsPoints,
        scanPoints,
        &deconvolvotronResult
        );
    QCOMPARE(e, eNoError);

    QCOMPARE(deconvolvotronResult.size(), pepStrVsPoints.size());

    //TODO fix pValue code

    QCOMPARE(deconvolvotronResult.at(0).first, "Pep1");
    QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.discScore, 0.5));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.pVal, 2.06302e-05));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(0).second.pValFrameFtest, 0.000492698));

    QCOMPARE(deconvolvotronResult.at(1).first, "Pep2");
    QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.discScore, 0.5));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.pVal, 3.78698e-05));
    // QVERIFY(MathUtils::tSame(deconvolvotronResult.at(1).second.pValFrameFtest, 0.000492698));

}


QTEST_MAIN(DeconvolvotronTests)

#include "DeconvolvotronTests.moc"
