#include "QuanTransitionRefinertron.h"

#include <QElapsedTimer>
#include <QtTest>

class QuanTransitionRefinertronTests : public QObject
{
    Q_OBJECT

public:
    QuanTransitionRefinertronTests() = default;
    ~QuanTransitionRefinertronTests() override = default;

private slots:

    static void smoothIntensities1Test();


};

void QuanTransitionRefinertronTests::smoothIntensities1Test() {

    // QSKIP("skipping for now");

    ERR_INIT

    const QVector<float> mzVals = {
        666.666, 666.667, 666.668, 666.686, 666.687, 666.688, 777.0
    };

    QHash<MzHashed, int> mzHashedVsCount;
    QuanTransitionRefinertron::groupTransitionsByMz(mzVals, 10.0, &mzHashedVsCount);

    QCOMPARE(mzHashedVsCount.value(666667), 3);
    QCOMPARE(mzHashedVsCount.value(666687), 3);
    QCOMPARE(mzHashedVsCount.value(777000), 1);

}


QTEST_MAIN(QuanTransitionRefinertronTests)

#include "QuanTransitionRefinertronTests.moc"
