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

    //TODO update this fragLib file
//    const QString &fragLibPath
//            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta.csv.fragLib");

    const double target = 454.957;
    const double targetWindowSize = 5.5;

    QMap<PeptideStringWithMods, MS2IonsSeparated> fragPreds;
    QMap<PeptideStringWithMods, bool> fragPredsIsDecoy;
    QMap<PeptideStringWithMods, double> fragPredsMass;
    QMap<PeptideStringWithMods, IRT> fragPredsIRT;

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibPath,
            2,
            3,
            target - targetWindowSize,
            target + targetWindowSize,
            &fragPreds,
            &fragPredsIsDecoy,
            &fragPredsMass
    );
    QCOMPARE(e, eNoError);

    FragLibIonRTree fragLibIonRTree;
    e = fragLibIonRTree.init(fragPreds);
    QCOMPARE(e, eNoError);

}

void FragLibIonRTreeTests::extractPointsTest() {

    ERR_INIT

    const QString fragLibPath
            = QStringLiteral("/home/anichols/Downloads/2022-04-23-decoys-contam-UP000005640_9606.target.20210602.human_plasma.fasta.csv.025.fragLib");

    //TODO update this fragLib file
//    const QString &fragLibPath
//            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta.csv.fragLib");

    const double target = 454.957;
    const double targetWindowSize = 5.5;

    QMap<PeptideStringWithMods, MS2IonsSeparated> fragPreds;
    QMap<PeptideStringWithMods, bool> fragPredsIsDecoy;
    QMap<PeptideStringWithMods, double> fragPredsMass;
    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibPath,
            2,
            3,
            target - targetWindowSize,
            target + targetWindowSize,
            &fragPreds,
            &fragPredsIsDecoy,
            &fragPredsMass
    );
    QCOMPARE(e, eNoError);

    FragLibIonRTree fragLibIonRTree;
    e = fragLibIonRTree.init(fragPreds);
    QCOMPARE(e, eNoError);

    const double ppmTolerance = 15.0;
    QMap<MzHashed, FrequencyPercent> mzHashVsFreqPct;
    e = fragLibIonRTree.addFrequencyPercentagesToFragLibIons(mzHashVsFreqPct);
    QCOMPARE(e, eNoError);
    QCOMPARE(static_cast<int>(mzHashVsFreqPct.value(147113) * 10000000), 97953);

}


QTEST_MAIN(FragLibIonRTreeTests)
#include "FragLibIonRTreeTests.moc"
