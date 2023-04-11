#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "AminoAcids.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PeptideMassRTree.h"


using namespace Error;


class PeptideMassRTreeTests : public QObject
{
    Q_OBJECT

public:

    PeptideMassRTreeTests() = default;
    ~PeptideMassRTreeTests() override = default;

private Q_SLOTS:

    void getPeptidesTest();

};

void PeptideMassRTreeTests::getPeptidesTest() {

    const QVector<PeptideStringWithMods> peps = {
        "VVVVVV",
        "V[10]VVVVV",
        "VVVTVV"
    };

    AminoAcids aminoAcids;

    PeptideMassRTree peptideMassRTree;
    Err e = peptideMassRTree.init(peps, aminoAcids);
    QCOMPARE(e, eNoError);

    QHash<PeptideStringWithMods , Mass> peptideStringWithModsTableVsMass;
    e = peptideMassRTree.getPeptides(
            610,
            614,
            &peptideStringWithModsTableVsMass
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideStringWithModsTableVsMass.size(), 1);
    QCOMPARE(peptideStringWithModsTableVsMass.contains(peps.at(0)), true);

    e = peptideMassRTree.getPeptides(
            620,
            624,
            &peptideStringWithModsTableVsMass
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideStringWithModsTableVsMass.size(), 1);
    QCOMPARE(peptideStringWithModsTableVsMass.contains(peps.at(1)), true);

    e = peptideMassRTree.getPeptides(
            615,
            614,
            &peptideStringWithModsTableVsMass
    );
    QCOMPARE(e, eValueError);

}


QTEST_MAIN(PeptideMassRTreeTests)

#include "PeptideMassRTreeTests.moc"
