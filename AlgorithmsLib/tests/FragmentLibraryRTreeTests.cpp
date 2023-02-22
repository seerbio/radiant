#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "FragmentLibraryRTree.h"


using namespace Error;


class FragmentLibraryRTreeTests : public QObject
{
    Q_OBJECT

public:

    FragmentLibraryRTreeTests() = default;
    ~FragmentLibraryRTreeTests() override = default;

private Q_SLOTS:

    void getFragmentIonsTest();

};

void FragmentLibraryRTreeTests::getFragmentIonsTest() {

    FragLibIon fli1;
    fli1.mzFrag = 1000.0;
    fli1.peptideMass = 998.993;
    fli1.peptideId = 666;

    FragLibIon fli2;
    fli2.mzFrag = 1000.01;
    fli2.peptideMass = 1000.0;
    fli2.peptideId = 6;

    FragLibIon fli3;
    fli3.mzFrag = 1000.0;
    fli3.peptideMass = 1001.0;
    fli3.peptideId = 66;

    FragLibIon fli4;
    fli4.mzFrag = 1000.0;
    fli4.peptideMass = 1997.99;
    fli4.peptideId = 6666;

    FragLibIon fli5;
    fli5.mzFrag = 1000.0;
    fli5.peptideMass = 2996.98;
    fli5.peptideId = 66666;

    FragLibIon fli6;
    fli6.mzFrag = 1000.011;
    fli6.peptideMass = 998.993;
    fli6.peptideId = 33;

    const QVector<FragLibIon> fragLibIons = {fli1, fli2, fli3, fli4, fli5, fli6, fli1};
    const double ppmTol = 10.0;
    const double targetWindow = 1.0;

    FragmentLibraryRTree fragmentLibraryRTree;
    Err e = fragmentLibraryRTree.init(fragLibIons, {1, 3}, ppmTol, targetWindow);
    QCOMPARE(e, eNoError);

    const QHash<PeptideId, double> t1 = fragmentLibraryRTree.getPeptidesTableIds(1000.0, 1000.0, {1.0,1.0});
    QCOMPARE(t1.size(), 4);

    const QVector<PeptideId> expectedKeys = {6666, 666, 6, 66666};
    for (PeptideId pi : expectedKeys) {
        const bool keyFound = t1.contains(pi);
        QCOMPARE(keyFound, true);
    }

}


QTEST_MAIN(FragmentLibraryRTreeTests)

#include "FragmentLibraryRTreeTests.moc"
