//
// Created by anichols on 11/07/2021.
//

#include "TandemSpectraDeconvolvotron.h"
#include "FragLibReader.h"

#include "ParallelUtils.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>

class TandemSpectraDeconvolvotronTests : public QObject
{

Q_OBJECT

public:
    TandemSpectraDeconvolvotronTests();
    ~TandemSpectraDeconvolvotronTests() override = default;


private Q_SLOTS:

    void deconvolveTandemSpectraTest();

};

TandemSpectraDeconvolvotronTests::TandemSpectraDeconvolvotronTests() : QObject(){}

void TandemSpectraDeconvolvotronTests::deconvolveTandemSpectraTest() {

    const MS2Ion i1(100.0, 1);
    const MS2Ion i2(200.0, 1);
    const MS2Ion i3(300.0, 1);

    const MS2Ion j1(100.0, 1);
    const MS2Ion j2(200.0, 1);
    const MS2Ion j3(600.0, 1);

    const QVector<QPointF> points = {
            {100, 1.5},
            {200, 1.5},
            {300, 1},
            {600, 0.5},
    };

    const QVector<MS2Ion> cand1 = {
            i1, i2, i3
    };

    const QVector<MS2Ion> cand2 = {
            j1, j2, j3
    };

    QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> input = {
            {"C1", cand1},
            {"C2", cand2}
    };

    QMap<PeptideSequenceChargeKey, double> result;

    TandemSpectraDeconvolvotron deconvolvotron;
    deconvolvotron.deconvolveTandemSpectra(
            points,
            input,
            &result
            );

    const QStringList expectedKeys = {"C1", "C2"};
    const QList<double> expectedVals = {1, 0.5};

    QCOMPARE(result.keys(), expectedKeys);
    QCOMPARE(result.value("C1"), 1);
    QCOMPARE(result.value("C2"), 0.5);

}

QTEST_MAIN(TandemSpectraDeconvolvotronTests)
#include "TandemSpectraDeconvolvotronTests.moc"
