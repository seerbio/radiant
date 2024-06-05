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
        {400.01, 10}
    };

    QMap<PeptideString, QVector<QPointF>> pepStrVsPoints = {
        {"Pep1", peptide1},
        {"Pep2", peptide2}
    };

    Deconvolvotron deconvolvotron;
    e = deconvolvotron.init(1);
    QCOMPARE(e, eNoError);

    QVector<QPair<IdStr, Score>> idStrVsScore;
    e = deconvolvotron.deconvolve(
        pepStrVsPoints,
        scanPoints,
        &idStrVsScore
        );
    QCOMPARE(e, eNoError);





}


QTEST_MAIN(DeconvolvotronTests)

#include "DeconvolvotronTests.moc"
