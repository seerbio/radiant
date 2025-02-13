#include <QtTest>
#include <QCoreApplication>

#include "GeneticFeatureSelection.h"



#include <iostream>

class GeneticFeatureSelectionTests : public QObject {
Q_OBJECT

public:
    GeneticFeatureSelectionTests() = default;

    ~GeneticFeatureSelectionTests() override = default;

private slots:

    static void OptimizorMsCalibratomaticOptimize();

};


void GeneticFeatureSelectionTests::OptimizorMsCalibratomaticOptimize() {
}

QTEST_MAIN(GeneticFeatureSelectionTests)

#include "GeneticFeatureSelectionTests.moc"
