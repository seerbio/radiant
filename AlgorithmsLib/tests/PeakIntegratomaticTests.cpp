#include "PeakIntegratomatic.h"

#include <QtTest>


class PeakIntegratomaticTests : public QObject
{
    Q_OBJECT

public:
    PeakIntegratomaticTests() = default;
    ~PeakIntegratomaticTests() override = default;

private slots:

    void initTest();

    void findAllPeakLimitsInXICTest();

    void troubleshoot();
};


void PeakIntegratomaticTests::initTest()
{
    ERR_INIT

    PeakIntegratomaticParameters params;
    params.filterLength = 1;

    PeakIntegratomatic peakIntegratomatic;

    e = peakIntegratomatic.init(params);
    QCOMPARE(e, eValueError);

    params.filterLength = 11;

    e = peakIntegratomatic.init(params);
    QCOMPARE(e, eNoError);

}


void PeakIntegratomaticTests::findAllPeakLimitsInXICTest() {

    ERR_INIT

    const QVector<double> xic = {
            198297, 16954, 21524, 17751, 13439, 12470, 41886, 54724, 76529, 53880, 50884, 14763, 72873,
            12431, 30320, 47730, 16619, 41425, 17111, 13134, 23280
    };

    PeakIntegratomaticParameters params;
    params.filterLength = 3;
    params.sigma = 1.0;
    params.signalToNoiseRatio = 1;
    params.smoothCount = 3;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(params);
    QCOMPARE(e, eNoError);

    QVector<PeakIntegrationIndexes> peakLimits;
    QVector<double> intensityVecSmoothed;
    e = peakIntegratomatic.findAllPeaksLimitsInXIC(
            xic,
            &peakLimits,
            &intensityVecSmoothed
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(peakLimits.size(), 1);
    QCOMPARE(peakLimits.back().first, 4);
    QCOMPARE(peakLimits.back().second, 23);

    e = peakIntegratomatic.findAllPeaksLimitsInXIC(
            {xic.front()},
            &peakLimits,
            &intensityVecSmoothed
    );
    QCOMPARE(e, eNoError);

    QCOMPARE(peakLimits.size(), 1);
    QCOMPARE(peakLimits.back().first, -3);
    QCOMPARE(peakLimits.back().second, 3);

}

void PeakIntegratomaticTests::troubleshoot() {

    ERR_INIT

    const QVector<double> xic = {
            683100, 637561, 1.94821e+06, 2.14361e+07, 6.92442e+07, 2.29084e+08, 4.92728e+08, 6.73543e+08, 8.58281e+08, 1.11204e+09, 1.02996e+09, 1.19564e+09, 1.41433e+09, 1.76781e+09, 1.42912e+09, 1.47188e+09, 1.28485e+09, 1.09102e+09, 9.37016e+08, 9.26473e+08, 6.27526e+08, 4.0335e+08, 2.47227e+08, 1.24627e+08, 6.19121e+07, 4.21142e+07, 3.58016e+07, 2.6134e+07, 2.73698e+07, 1.8055e+07, 1.30663e+07, 1.19533e+07, 6.38015e+06, 9.94669e+06, 8.31145e+06, 6.02947e+06, 8.96082e+06, 5.35547e+06, 8.73866e+06, 6.71771e+06, 1.04763e+07, 6.63395e+06, 6.11291e+06, 4.67784e+06, 6.85373e+06, 3.71956e+06, 6.09384e+06, 4.9013e+06, 4.37596e+06, 4.80583e+06, 2.65859e+06, 4.61011e+06, 5.29583e+06, 4.58558e+06, 4.5843e+06, 3.77166e+06, 4.0135e+06, 4.7433e+06, 5.45258e+06, 4.75796e+06, 4.23713e+06, 2.42815e+06, 2.88705e+06, 3.26178e+06, 5.24376e+06, 1.34047e+06, 3.24541e+06, 2.23853e+06, 2.93705e+06, 2.08196e+06, 3.51306e+06, 3.07677e+06, 2.04034e+06, 2.1163e+06, 2.27004e+06, 1.68733e+06, 2.32565e+06, 2.32943e+06, 1.01183e+06, 2.37571e+06, 1.62404e+06, 1.2743e+06, 2.3072e+06, 2.57018e+06, 3.29065e+06, 1.84329e+06, 2.36236e+06, 1.76459e+06, 1.32381e+06, 1.23579e+06, 2.60593e+06, 2.42852e+06, 2.18522e+06, 2.85201e+06, 2.57898e+06, 2.65825e+06, 1.88165e+06, 1.18385e+06, 1.04629e+06, 1.95062e+06, 1.35698e+06, 836220, 693558, 2.00161e+06, 2.03359e+06, 1.75181e+06, 973249, 1.0736e+06, 1.64345e+06, 792523, 1.35283e+06, 1.48632e+06, 2.02371e+06, 1.3347e+06, 1.97834e+06, 815144, 1.78048e+06, 644160, 954895, 945109, 1.49231e+06, 1.40992e+06, 1.32281e+06, 783219, 704865, 3.19979e+06, 796806, 1.63386e+06, 1.01376e+06, 1.97446e+06, 1.85681e+06, 2.29583e+06, 1.37941e+06, 2.17791e+06, 1.64737e+06, 1.03536e+06, 1.10366e+06, 537947, 851528, 1.58797e+06, 1.7302e+06, 1.17549e+06, 535842, 781794, 752733, 764563, 1.44993e+06, 1.24531e+06, 421828, 659686, 1.25447e+06, 779897, 1.10535e+06, 712997, 858159, 479404, 1.09687e+06, 1.04983e+06, 1.39503e+06, 558596, 1.21065e+06, 775444, 957731, 666247, 637820, 983087, 1.25676e+06, 335943, 771547, 526644, 758424, 1.02261e+06, 550738, 1.2414e+06, 407835, 275116, 973569, 704034, 317924, 1.13804e+06, 871434, 596976, 619573, 588770, 1.0324e+06, 665885, 592333, 1.06875e+06, 305434, 679738, 460539, 542297, 1.3145e+06, 901281, 899041, 897429, 945524, 824468, 1.21035e+06, 553438, 1.05773e+06, 627385, 604191, 418294, 578730, 532376, 532095, 574255, 197265, 406817, 863067, 582606, 339582, 557238, 603349, 468139, 674008, 448096, 645076, 665035, 798168, 696605, 676162, 475404, 326880, 537244, 441279, 888373, 368966, 469136, 444062, 502362, 583502, 634340, 313966, 449843, 363408, 378235, 611000, 678750, 691138, 535717, 599261, 351143, 456993, 339499, 735600, 584803, 538955, 838977, 478828, 539276, 517798, 375168, 687691, 458030, 371398, 520492, 674733, 925356, 650800, 789707, 723242, 624297, 308958, 222706, 263559, 683802, 493047, 1.30123e+06, 774124, 1.21462e+06, 1.04835e+06, 558385, 1.31298e+06, 691405, 462973, 635231, 681066, 752560, 641209, 718388, 648227, 552868, 814124, 2.028e+06, 1.66319e+06
    };

    PeakIntegratomaticParameters params;
    params.filterLength = 3;
    params.sigma = 1.0;
    params.signalToNoiseRatio = 1;
    params.smoothCount = 3;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(params);
    QCOMPARE(e, eNoError);

    QVector<PeakIntegrationIndexes> peakLimits;
    QVector<double> intensityVecSmoothed;
    e = peakIntegratomatic.findAllPeaksLimitsInXIC(
            xic,
            &peakLimits,
            &intensityVecSmoothed
    );
    QCOMPARE(e, eNoError);
}

QTEST_MAIN(PeakIntegratomaticTests)

#include "PeakIntegratomaticTests.moc"
