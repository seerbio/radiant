//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "FragLibIonRTree.h"
#include "FragLibReader.h"

#include <QtTest/QtTest>

class FragLibIonRTreeTests : public QObject
{
    Q_OBJECT

public:
    FragLibIonRTreeTests() = default;
    ~FragLibIonRTreeTests() override = default;

private Q_SLOTS:

    void initTest();
    void extractPointsTest();

};


void FragLibIonRTreeTests::initTest() {

    ERR_INIT

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Downloads/2022-04-23-decoys-contam-UP000005640_9606.target.20210602.human_plasma.fasta.csv.025.fragLib");

    const double target = 454.957;
    const double targetWindowSize = 5.5;

    QMap<PeptideStringWithMods, MS2IonsSeparated> fragPreds;
    QMap<PeptideStringWithMods, bool> fragPredsIsDecoy;

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibPath,
            2,
            3,
            target - targetWindowSize,
            target + targetWindowSize,
            &fragPreds,
            &fragPredsIsDecoy
    );
    QCOMPARE(e, eNoError);

    FragLibIonRTree fragLibIonRTree;
    e = fragLibIonRTree.init(fragPreds);
    QCOMPARE(e, eNoError);

}


void FragLibIonRTreeTests::extractPointsTest() {

    ERR_INIT



}


QTEST_MAIN(FragLibIonRTreeTests)
#include "FragLibIonRTreeTests.moc"
