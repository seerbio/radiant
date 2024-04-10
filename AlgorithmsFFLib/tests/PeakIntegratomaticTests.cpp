#include "PeakIntegratomatic.h"

#include <QHash>
#include <QtTest>

class PeakIntegratomaticTests : public QObject
{
    Q_OBJECT

public:
    PeakIntegratomaticTests() = default;
    ~PeakIntegratomaticTests() override = default;

private slots:

    static void buildPeakIntegratomaticParamsTest();

    static void initTest();
    static void simpleIntegratorTest();

};


void PeakIntegratomaticTests::buildPeakIntegratomaticParamsTest() {

    PeakIntegratomaticParameters params = PeakIntegratomaticParameters::buildPeakIntegratomaticParams(
            PythiaParameterReader::genericPythiaParametersForTests()
            );

    QCOMPARE(params.filterLength, 5);
    QCOMPARE(params.smoothCount, 2);
    QCOMPARE(params.sigma, 1.0);
    QCOMPARE(params.signalToNoiseRatio, 2.0);
    QCOMPARE(MathUtils::pRound(static_cast<double>(params.stopThresholdFraction), 1), static_cast<double>(0.2));

}

void PeakIntegratomaticTests::initTest() {

    ERR_INIT

    PeakIntegratomaticParameters paramsBad;
    paramsBad.filterLength = -1;

    PeakIntegratomatic peakIntegratomatic;

    e = peakIntegratomatic.init(paramsBad);
    QCOMPARE(e, eValueError);

    PeakIntegratomaticParameters paramsGood = PeakIntegratomaticParameters::buildPeakIntegratomaticParams(
            PythiaParameterReader::genericPythiaParametersForTests()
    );

    e = peakIntegratomatic.init(paramsGood);
    QCOMPARE(e, eNoError);

}

void PeakIntegratomaticTests::simpleIntegratorTest() {

    ERR_INIT

    const QVector<float> testVec = {
            1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
            1, 2, 2, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 2, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 3, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 1, 0, 0, 0, 0, 1, 0, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2,
            4, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 1, 0, 1, 2, 2, 2, 0, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 1, 0, 0, 0, 0, 1, 0, 3, 1, 1, 2, 4, 2, 0, 0, 0, 1, 2, 0, 2, 2,
            3, 1, 0, 1, 0, 1, 0, 3, 4, 5, 10, 12, 12, 11, 10, 10, 7, 5, 2, 3, 2, 1, 1, 2, 2, 2, 3, 2, 4, 11, 12,
            12, 12, 12, 12, 11, 11, 11, 11, 11, 10, 10, 6, 3, 5, 4, 4, 3, 5, 4, 3, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1,
            3, 1, 1, 2, 1, 1, 2, 4, 2, 1, 2, 0, 0, 1, 0, 1, 2, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0
    };

    PeakIntegratomaticParameters peakIntegratomaticParameters;
    peakIntegratomaticParameters.filterLength = 3;
    peakIntegratomaticParameters.smoothCount = 1;
    peakIntegratomaticParameters.sigma = 1.0;
    peakIntegratomaticParameters.signalToNoiseRatio = 2.0;
    peakIntegratomaticParameters.stopThresholdFraction = 0.2;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(peakIntegratomaticParameters);
    QCOMPARE(e, eNoError);


    QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexes;
    e = peakIntegratomatic.simpleIntegrator(
            testVec,
            10,
            5,
            &peakIntegrationIndexes
            );
    QCOMPARE(e, eNoError);

    using PII = QPair<PeakIntegrationIndexes, float>;

    std::sort(
            peakIntegrationIndexes.rbegin(),
            peakIntegrationIndexes.rend(),
            [](const PII &l, const PII &r){return l.second <r.second;}
            );

    QCOMPARE(peakIntegrationIndexes.size(), 5);

    QCOMPARE(peakIntegrationIndexes.front().first, QPair(256, 274));
    QCOMPARE( static_cast<int>(std::round(peakIntegrationIndexes.front().second)), 21);

    QCOMPARE(peakIntegrationIndexes.at(1).first, QPair(237, 251));
    QCOMPARE( static_cast<int>(std::round(peakIntegrationIndexes.at(1).second)), 20);

    QCOMPARE(peakIntegrationIndexes.back().first, QPair(297, 301));
    QCOMPARE( static_cast<int>(std::round(peakIntegrationIndexes.back().second)), 5);

}


QTEST_MAIN(PeakIntegratomaticTests)

#include "PeakIntegratomaticTests.moc"
