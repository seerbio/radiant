#include <QtTest>
#include <QCoreApplication>

#include "IsotopicDistributionBuilder.h"

#include <QVector>

class IsotopeDistributionBuilderTests : public QObject
{
    Q_OBJECT

public:
    IsotopeDistributionBuilderTests();
    ~IsotopeDistributionBuilderTests();

private slots:
    void getIsotopicDistributionTest1();
    void getIsotopicDistributionTest2();

};

IsotopeDistributionBuilderTests::IsotopeDistributionBuilderTests(){}

IsotopeDistributionBuilderTests::~IsotopeDistributionBuilderTests(){}

void IsotopeDistributionBuilderTests::getIsotopicDistributionTest1()
{


    const QVector<double> dist
            = IsotopicDistributionBuilder::getIsotopicDistribution(QStringLiteral("QVTLRESGPALVKPTQT").split(""));

    QVector<double> expected
            = {1, 0.986017, 0.507067, 0.180479, 0.0498093, 0.0113286, 0.0022049};

    QCOMPARE(dist.size(), expected.size());

    const int muliplier = 100;
    for(int i = 0; i < dist.size(); ++i){
        QCOMPARE(std::round(dist.at(i) * muliplier), std::round(expected.at(i) * muliplier));
    }

    const QVector<double> distCys
            = IsotopicDistributionBuilder::getIsotopicDistribution(QStringLiteral("CCCC").split(""));

}


void IsotopeDistributionBuilderTests::getIsotopicDistributionTest2()
{

    const QVector<double> dist = IsotopicDistributionBuilder::getIsotopicDistribution(1806.0);

    QVector<double> expected
        ={1, 0.991266, 0.549142, 0.219717, 0.0693163, 0.0180092, 0.00397333, 0.000764337};

    QCOMPARE(dist.size(), expected.size());

    const int muliplier = 100;
    for(int i = 0; i < dist.size(); ++i){
        QCOMPARE(std::round(dist.at(i) * muliplier), std::round(expected.at(i) * muliplier));
    }

}


QTEST_MAIN(IsotopeDistributionBuilderTests)

#include "IsotopeDistributionBuilderTests.moc"
