#include "BiophysicalCalcs.h"

#include "MS2Ion.h"

#include <QHash>
#include <QtTest>

class MS2IonTests : public QObject
{
    Q_OBJECT

public:
    MS2IonTests() = default;
    ~MS2IonTests() override = default;

private slots:

    static void getIonLabelInfoTest();

    static void sortMS2IonsMzAscTest();

    static void sortMS2IonsIntensityDescTest();

    static void filterMS2IonsByMzTest();

    static void ms2IonsToScanPointsTests();

private:

    static QVector<MS2Ion> buildMS2Ions();

};

QVector<MS2Ion> MS2IonTests::buildMS2Ions() {

    MS2Ion ms2Ion1;
    ms2Ion1.mz = 888.8;
    ms2Ion1.intensity = 1.0;

    MS2Ion ms2Ion2;
    ms2Ion2.mz = 777.7;
    ms2Ion2.intensity = 0.5;

    MS2Ion ms2Ion3;
    ms2Ion3.mz = 666.6;
    ms2Ion3.intensity = 0.75;

    MS2Ion ms2Ion4;
    ms2Ion4.mz = 999.9;
    ms2Ion4.intensity = 0.25;

    QVector<MS2Ion> ms2Ions = {
            ms2Ion1,
            ms2Ion2,
            ms2Ion3,
            ms2Ion4
    };

    return ms2Ions;
}

void MS2IonTests::getIonLabelInfoTest() {

    MS2Ion ms2Ion;
    ms2Ion.mz = 666.6;
    ms2Ion.intensity = 1.0;
    ms2Ion.charge = 2;

    QPair<IonIndex, IonType> ionInfo;

    ms2Ion.ionLabel = QStringLiteral("b13^2");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "b^2");

    ms2Ion.ionLabel = QStringLiteral("b13");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "b");

    ms2Ion.ionLabel = QStringLiteral("y13^2");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "y^2");

    ms2Ion.ionLabel = QStringLiteral("y13");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "y");

    ms2Ion.ionLabel = QStringLiteral("y13-NH3");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "y-NH3");

    ms2Ion.ionLabel = QStringLiteral("b13-H2O");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "b-H2O");

    ms2Ion.ionLabel = QStringLiteral("a13");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 13);
    QCOMPARE(ionInfo.second, "a");

    ms2Ion.ionLabel = QStringLiteral("p-NH3");
    ms2Ion.getIonLabelInfo(&ionInfo);
    QCOMPARE(ionInfo.first, 0);
    QCOMPARE(ionInfo.second, "p-NH3");

}

void MS2IonTests::sortMS2IonsMzAscTest() {

    QVector<MS2Ion> ms2Ions = buildMS2Ions();

    MS2Ion::sortMS2IonsMzAsc(&ms2Ions);

    QVERIFY(MathUtils::tSame(ms2Ions.at(0).mz, 666.6f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(1).mz, 777.7f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(2).mz, 888.8f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(3).mz, 999.9f));

}

void MS2IonTests::sortMS2IonsIntensityDescTest() {

    QVector<MS2Ion> ms2Ions = buildMS2Ions();

    MS2Ion::sortMS2IonsIntensityDesc(&ms2Ions);

    QVERIFY(MathUtils::tSame(ms2Ions.at(0).mz, 888.8f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(1).mz, 666.6f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(2).mz, 777.7f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(3).mz, 999.9f));

}

void MS2IonTests::filterMS2IonsByMzTest() {

    QVector<MS2Ion> ms2Ions = buildMS2Ions();

    MS2Ion::filterMS2IonsByMz(700.0, 900.0, &ms2Ions);
    MS2Ion::sortMS2IonsMzAsc(&ms2Ions);

    QVERIFY(MathUtils::tSame(ms2Ions.at(0).mz, 777.7f));
    QVERIFY(MathUtils::tSame(ms2Ions.at(1).mz, 888.8f));

}

void MS2IonTests::ms2IonsToScanPointsTests() {

    const QVector<MS2Ion> ms2Ions = buildMS2Ions();
    ScanPoints scanPoints = MS2Ion::ms2IonsToScanPoints(ms2Ions);

    QVERIFY(MathUtils::tSame(scanPoints.at(0).x(), 888.8f));
    QVERIFY(MathUtils::tSame(scanPoints.at(1).x(), 777.7f));
    QVERIFY(MathUtils::tSame(scanPoints.at(2).x(), 666.6f));
    QVERIFY(MathUtils::tSame(scanPoints.at(3).x(), 999.9f));

}


QTEST_MAIN(MS2IonTests)

#include "MS2IonTests.moc"
