#include "ConvexOptimizer.h"

#include <QElapsedTimer>
#include <QtTest/QtTest>


class ConvexOptimizerTests : public QObject
{
    Q_OBJECT

public:
    ConvexOptimizerTests() = default;
    ~ConvexOptimizerTests() = default;

private Q_SLOTS:

    void loadModelTest();


};

void ConvexOptimizerTests::loadModelTest() {

    ERR_INIT

    ConvexOptimizer convexOptimizer;
    e = convexOptimizer.test();

}


QTEST_MAIN(ConvexOptimizerTests)
#include "ConvexOptimizerTests.moc"
