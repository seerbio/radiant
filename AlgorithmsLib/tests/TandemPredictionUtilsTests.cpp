#include "BiophysicalCalcs.h"

#include "AminoAcids.h"
#include "Error.h"
#include "TandemPredictionUtils.h"

#include <QtTest>

class TandemPredictionUtilsTests : public QObject
{
    Q_OBJECT

public:
    TandemPredictionUtilsTests() = default;
    ~TandemPredictionUtilsTests() override = default;

private slots:

    void buildTandemFragmentMassesTest();
//    void test();

};


void TandemPredictionUtilsTests::buildTandemFragmentMassesTest()
{
    const int charge = 2;
    const QStringList ionListCharge2 = TandemPredictionUtils::buildIonLabels(charge);

    const QStringList ionListCharge2Expected = {"p", "p-H2O", "p-NH3", "y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8", "y9", "y10", "y11", "y12",
             "y13", "y14", "y15", "y16", "y17", "y18", "y19", "y20", "y21", "y22", "y23", "y24", "y25", "y26", "y27", "y28", "y29", "y30", "y31", "y32",
             "y33", "y34", "y35", "y36", "y37", "y38", "y39", "y40", "y1^2", "y2^2", "y3^2", "y4^2", "y5^2", "y6^2", "y7^2", "y8^2", "y9^2", "y10^2",
             "y11^2", "y12^2", "y13^2", "y14^2", "y15^2", "y16^2", "y17^2", "y18^2", "y19^2", "y20^2", "y21^2", "y22^2", "y23^2", "y24^2", "y25^2",
             "y26^2", "y27^2", "y28^2", "y29^2", "y30^2", "y31^2", "y32^2", "y33^2", "y34^2", "y35^2", "y36^2", "y37^2", "y38^2", "y39^2", "y40^2",
             "y1-NH3", "y2-NH3", "y3-NH3", "y4-NH3", "y5-NH3", "y6-NH3", "y7-NH3", "y8-NH3", "y9-NH3", "y10-NH3", "y11-NH3", "y12-NH3",
             "y13-NH3", "y14-NH3", "y15-NH3", "y16-NH3", "y17-NH3", "y18-NH3", "y19-NH3", "y20-NH3", "y21-NH3", "y22-NH3", "y23-NH3",
             "y24-NH3", "y25-NH3", "y26-NH3", "y27-NH3", "y28-NH3", "y29-NH3", "y30-NH3", "y31-NH3", "y32-NH3", "y33-NH3", "y34-NH3",
            "y35-NH3", "y36-NH3", "y37-NH3", "y38-NH3", "y39-NH3", "y40-NH3", "y1-H2O", "y2-H2O", "y3-H2O", "y4-H2O", "y5-H2O", "y6-H2O",
            "y7-H2O", "y8-H2O", "y9-H2O", "y10-H2O", "y11-H2O", "y12-H2O", "y13-H2O", "y14-H2O", "y15-H2O", "y16-H2O", "y17-H2O",
            "y18-H2O", "y19-H2O", "y20-H2O", "y21-H2O", "y22-H2O", "y23-H2O", "y24-H2O", "y25-H2O", "y26-H2O", "y27-H2O", "y28-H2O",
            "y29-H2O", "y30-H2O", "y31-H2O", "y32-H2O", "y33-H2O", "y34-H2O", "y35-H2O", "y36-H2O", "y37-H2O", "y38-H2O", "y39-H2O",
            "y40-H2O", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19",
            "b20", "b21", "b22", "b23", "b24", "b25", "b26", "b27", "b28", "b29", "b30", "b31", "b32", "b33", "b34", "b35", "b36", "b37", "b38",
            "b39", "b40", "b1^2", "b2^2", "b3^2", "b4^2", "b5^2", "b6^2", "b7^2", "b8^2", "b9^2", "b10^2", "b11^2", "b12^2", "b13^2", "b14^2",
            "b15^2", "b16^2", "b17^2", "b18^2", "b19^2", "b20^2", "b21^2", "b22^2", "b23^2", "b24^2", "b25^2", "b26^2", "b27^2", "b28^2", "b29^2",
            "b30^2", "b31^2", "b32^2", "b33^2", "b34^2", "b35^2", "b36^2", "b37^2", "b38^2", "b39^2", "b40^2", "b1-NH3", "b2-NH3", "b3-NH3",
            "b4-NH3", "b5-NH3", "b6-NH3", "b7-NH3", "b8-NH3", "b9-NH3", "b10-NH3", "b11-NH3", "b12-NH3", "b13-NH3", "b14-NH3", "b15-NH3",
            "b16-NH3", "b17-NH3", "b18-NH3", "b19-NH3", "b20-NH3", "b21-NH3", "b22-NH3", "b23-NH3", "b24-NH3", "b25-NH3", "b26-NH3",
            "b27-NH3", "b28-NH3", "b29-NH3", "b30-NH3", "b31-NH3", "b32-NH3", "b33-NH3", "b34-NH3", "b35-NH3", "b36-NH3", "b37-NH3",
            "b38-NH3", "b39-NH3", "b40-NH3", "b1-H2O", "b2-H2O", "b3-H2O", "b4-H2O", "b5-H2O", "b6-H2O", "b7-H2O", "b8-H2O", "b9-H2O",
            "b10-H2O", "b11-H2O", "b12-H2O", "b13-H2O", "b14-H2O", "b15-H2O", "b16-H2O", "b17-H2O", "b18-H2O", "b19-H2O", "b20-H2O", "b21-H2O",
            "b22-H2O", "b23-H2O", "b24-H2O", "b25-H2O", "b26-H2O", "b27-H2O", "b28-H2O", "b29-H2O", "b30-H2O", "b31-H2O", "b32-H2O",
            "b33-H2O", "b34-H2O", "b35-H2O", "b36-H2O", "b37-H2O", "b38-H2O", "b39-H2O", "b40-H2O", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8",
            "a9", "a10", "a11", "a12", "a13", "a14", "a15", "a16", "a17", "a18", "a19", "a20", "a21", "a22", "a23", "a24", "a25", "a26", "a27",
            "a28", "a29", "a30", "a31", "a32", "a33", "a34", "a35", "a36", "a37", "a38", "a39", "a40"};

    QCOMPARE(ionListCharge2.size(), ionListCharge2Expected.size());
    QCOMPARE(ionListCharge2, ionListCharge2Expected);

}


QTEST_MAIN(TandemPredictionUtilsTests)

#include "TandemPredictionUtilsTests.moc"
