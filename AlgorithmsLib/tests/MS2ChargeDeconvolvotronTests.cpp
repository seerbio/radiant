#include "BiophysicalCalcs.h"
#include "Error.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsUtils.h"
#include "ParallelUtils.h"

#include <QtTest>


class MS2ChargeDeconvolvotronTests : public QObject
{
    Q_OBJECT

public:
    MS2ChargeDeconvolvotronTests() = default;
    ~MS2ChargeDeconvolvotronTests() override = default;

private slots:

    void initTest();
    void buildMzValsForChargeMonoTest();
    void predictMS2ChargeAndMonoOffsetTest();
    void deisotopeScanPointsTest();
    void deisotopeScanPointsRealDataTest();


private:

    ScanPoints buildToyScanPoints();

private:

    const QVector<ScanPoint> m_scanPoints = {
            QPointF(145.061, 540509), QPointF(147.113, 32804872), QPointF(148.06, 517785), QPointF(148.116, 2306725), QPointF(149.045, 508699), QPointF(150.056, 1046070),
            QPointF(152.107, 420139), QPointF(153.068, 281177), QPointF(155.118, 6388912), QPointF(156.121, 471498), QPointF(157.097, 1246425), QPointF(158.081, 19641120),
            QPointF(159.085, 1364219), QPointF(161.093, 4479668), QPointF(162.095, 380779), QPointF(167.082, 10758629), QPointF(167.091, 230739), QPointF(168.065, 5092302),
            QPointF(168.085, 405403), QPointF(169.061, 4061765), QPointF(169.069, 452704), QPointF(169.098, 608747), QPointF(169.133, 602641), QPointF(171.077, 3074904),
            QPointF(173.128, 63788924), QPointF(174.132, 5790877), QPointF(175.119, 392805), QPointF(175.135, 301513), QPointF(181.062, 396186), QPointF(183.113, 2492807),
            QPointF(184.107, 1597685), QPointF(185.092, 62008724), QPointF(186.077, 395667), QPointF(186.095, 5628402), QPointF(186.124, 591486), QPointF(187.097, 590183),
            QPointF(187.144, 2016346), QPointF(188.147, 281227), QPointF(189.086, 5124722), QPointF(190.091, 396439), QPointF(195.077, 18421740), QPointF(196.08, 2161968),
            QPointF(197.092, 478277), QPointF(197.128, 3479838), QPointF(198.123, 1812320), QPointF(198.134, 265367), QPointF(199.108, 445803), QPointF(199.127, 323378),
            QPointF(201.123, 29620824), QPointF(202.127, 2959791), QPointF(203.103, 36450804), QPointF(204.105, 3106530), QPointF(205.106, 322162), QPointF(208.11, 841872),
            QPointF(209.092, 5423484), QPointF(209.104, 420520), QPointF(209.141, 1494142), QPointF(210.097, 321894), QPointF(213.087, 32766430), QPointF(214.091, 3750324),
            QPointF(214.119, 4616792), QPointF(215.102, 494838), QPointF(215.125, 435320), QPointF(215.138, 6801422), QPointF(216.132, 346508), QPointF(216.144, 445692),
            QPointF(217.081, 399585), QPointF(218.15, 630735), QPointF(219.134, 317983), QPointF(221.056, 374991), QPointF(222.133, 711361), QPointF(225.134, 953825),
            QPointF(226.118, 24889962), QPointF(226.645, 1586472), QPointF(226.66, 2736608), QPointF(227.103, 1688874), QPointF(227.122, 2338423), QPointF(228.133, 541013),
            QPointF(231.081, 427652), QPointF(231.096, 13043223), QPointF(231.113, 364686), QPointF(231.138, 552758), QPointF(232.101, 1489286), QPointF(232.128, 358724),
            QPointF(233.131, 700499), QPointF(236.103, 581788), QPointF(237.086, 608571), QPointF(237.134, 825019), QPointF(238.081, 2068305), QPointF(238.154, 607729),
            QPointF(239.065, 374198), QPointF(239.139, 534998), QPointF(243.145, 9997223), QPointF(244.111, 340705), QPointF(244.129, 27387556), QPointF(244.145, 816363),
            QPointF(245.133, 3035518), QPointF(246.131, 407340), QPointF(246.18, 352528), QPointF(252.134, 302350), QPointF(254.099, 192490), QPointF(254.114, 11343426),
            QPointF(254.131, 488207), QPointF(255.099, 554277), QPointF(255.116, 1358920), QPointF(256.092, 2684250), QPointF(256.166, 1007779), QPointF(257.096, 299426),
            QPointF(257.15, 418744), QPointF(258.143, 517250), QPointF(261.14, 2320448), QPointF(261.155, 124303496), QPointF(261.191, 829321), QPointF(262.159, 13853521),
            QPointF(263.161, 1239562), QPointF(264.098, 7157454), QPointF(265.101, 1081510), QPointF(266.148, 984171), QPointF(267.133, 390786), QPointF(267.151, 395205),
            QPointF(268.128, 271161), QPointF(268.202, 455511), QPointF(270.145, 616731), QPointF(270.161, 269351), QPointF(272.122, 2905467), QPointF(273.127, 329747),
            QPointF(274.176, 716403), QPointF(278.674, 2757820), QPointF(279.177, 931294), QPointF(279.649, 329946), QPointF(280.129, 295517), QPointF(281.125, 1360513),
            QPointF(282.109, 12346555), QPointF(283.112, 1783748), QPointF(284.16, 17254720), QPointF(285.163, 1982076), QPointF(287.171, 591667), QPointF(287.682, 5154351),
            QPointF(288.184, 1528528), QPointF(288.683, 382139), QPointF(294.179, 321369), QPointF(296.197, 2443925), QPointF(297.154, 511797), QPointF(297.2, 277940),
            QPointF(298.14, 337435), QPointF(299.134, 515816), QPointF(300.117, 14491142), QPointF(301.122, 2026043), QPointF(302.168, 13946829), QPointF(303.173, 2272265),
            QPointF(304.177, 259616), QPointF(304.671, 930774), QPointF(313.16, 418586), QPointF(313.186, 1130880), QPointF(313.678, 1394036), QPointF(314.18, 527011),
            QPointF(314.207, 411023), QPointF(315.141, 361457), QPointF(315.165, 9834561), QPointF(316.151, 3514281), QPointF(316.171, 1147116), QPointF(317.152, 680242),
            QPointF(318.13, 1920668), QPointF(320.195, 308837), QPointF(322.192, 22121636), QPointF(322.22, 370278), QPointF(322.696, 6696398), QPointF(323.173, 2907844),
            QPointF(323.197, 1649254), QPointF(324.176, 358565), QPointF(324.201, 281867), QPointF(325.187, 2995150), QPointF(326.193, 479520), QPointF(331.123, 326517),
            QPointF(331.199, 1771360), QPointF(331.699, 557394), QPointF(333.178, 607417), QPointF(338.205, 307284), QPointF(341.018, 554899), QPointF(341.177, 423845),
            QPointF(342.211, 1032826), QPointF(343.196, 4756052), QPointF(343.22, 266087), QPointF(344.2, 895159), QPointF(344.372, 272039), QPointF(349.185, 13067234),
            QPointF(350.167, 244107), QPointF(350.192, 2390394), QPointF(350.702, 2084442), QPointF(351.203, 370968), QPointF(353.184, 842301), QPointF(354.706, 819190),
            QPointF(355.208, 389539), QPointF(355.236, 3009764), QPointF(355.699, 424190), QPointF(356.239, 336313), QPointF(357.218, 385958), QPointF(359.027, 794397),
            QPointF(360.226, 48515232), QPointF(361.226, 8902792), QPointF(362.229, 1257608), QPointF(363.71, 5226366), QPointF(364.212, 2230376), QPointF(364.704, 572733),
            QPointF(365.218, 376769), QPointF(367.199, 9842840), QPointF(368.203, 1947598), QPointF(370.207, 540077), QPointF(371.192, 1058218), QPointF(372.714, 2018054),
            QPointF(373.221, 478950), QPointF(373.243, 1237404), QPointF(377.181, 2099837), QPointF(378.181, 302893), QPointF(381.719, 1101588), QPointF(382.221, 463602),
            QPointF(383.227, 9557117), QPointF(384.233, 2264440), QPointF(385.207, 16964968), QPointF(386.211, 3528108), QPointF(387.211, 500334), QPointF(392.703, 419010),
            QPointF(395.194, 6895969), QPointF(396.195, 1560508), QPointF(397.216, 744271), QPointF(398.206, 716668), QPointF(401.241, 1955195), QPointF(401.709, 798837),
            QPointF(402.207, 831499), QPointF(402.241, 546082), QPointF(403.22, 1555303), QPointF(406.222, 1583197), QPointF(406.719, 523150), QPointF(408.225, 746104),
            QPointF(410.22, 2274034), QPointF(410.72, 1441491), QPointF(411.221, 875475), QPointF(413.205, 7054888), QPointF(414.21, 1275117), QPointF(414.236, 483247),
            QPointF(417.248, 1956185), QPointF(418.252, 473514), QPointF(419.228, 3587473), QPointF(419.728, 1822552), QPointF(420.23, 430902), QPointF(425.251, 454458),
            QPointF(426.236, 1576046), QPointF(427.237, 389538), QPointF(428.231, 4110869), QPointF(428.731, 1810174), QPointF(429.199, 301094), QPointF(429.234, 1643585),
            QPointF(434.204, 1547391), QPointF(441.251, 648211), QPointF(443.262, 4688828), QPointF(443.754, 350375), QPointF(444.246, 5541368), QPointF(445.25, 1692983),
            QPointF(446.262, 570400), QPointF(446.295, 370670), QPointF(447.202, 1937517), QPointF(448.254, 375562), QPointF(450.237, 814139), QPointF(452.215, 951218),
            QPointF(452.25, 358665), QPointF(453.244, 283245), QPointF(460.222, 2509588), QPointF(461.189, 862181), QPointF(461.237, 2195084), QPointF(461.272, 166405200),
            QPointF(462.273, 37995240), QPointF(463.277, 6424731), QPointF(464.211, 344057), QPointF(464.283, 633416), QPointF(466.265, 1070324), QPointF(468.247, 1872172),
            QPointF(469.244, 382568), QPointF(470.221, 1068057), QPointF(470.258, 325177), QPointF(471.259, 1796534), QPointF(472.26, 428350), QPointF(478.229, 5983044),
            QPointF(479.226, 1470384), QPointF(484.273, 896528), QPointF(485.279, 379022), QPointF(486.256, 533080), QPointF(496.241, 22249630), QPointF(497.239, 4747440),
            QPointF(498.247, 1564937), QPointF(502.294, 411408), QPointF(514.252, 15596475), QPointF(515.248, 2989522), QPointF(515.36, 299023), QPointF(516.26, 802853),
            QPointF(516.295, 464541), QPointF(517.267, 306243), QPointF(530.332, 305720), QPointF(531.294, 498495), QPointF(539.318, 480766), QPointF(540.274, 622612),
            QPointF(549.3, 802941), QPointF(551.276, 386282), QPointF(556.345, 5013958), QPointF(556.401, 316647), QPointF(557.343, 1596685), QPointF(558.344, 590148),
            QPointF(559.287, 1030687), QPointF(567.311, 2946755), QPointF(568.31, 781256), QPointF(574.357, 40627264), QPointF(575.358, 12658206), QPointF(576.361, 2486624),
            QPointF(577.302, 437217), QPointF(584.342, 617850), QPointF(585.335, 616767), QPointF(595.313, 2256378), QPointF(596.316, 849868), QPointF(608.343, 903826),
            QPointF(613.315, 756789), QPointF(617.364, 318702), QPointF(623.27, 394700), QPointF(625.358, 467997), QPointF(626.269, 723655), QPointF(626.345, 1120085),
            QPointF(640.318, 525136), QPointF(643.372, 10090623), QPointF(644.375, 4609924), QPointF(645.364, 1622897), QPointF(646.36, 295878), QPointF(649.37, 884580),
            QPointF(650.364, 294489), QPointF(661.239, 698940), QPointF(661.381, 84275312), QPointF(662.393, 29054324), QPointF(663.386, 7162041), QPointF(664.393, 724569),
            QPointF(671.376, 407189), QPointF(709.35, 460774), QPointF(710.338, 529019), QPointF(726.407, 694137), QPointF(727.407, 500607), QPointF(744.426, 1267274),
            QPointF(745.421, 563406), QPointF(748.433, 631572), QPointF(754.403, 604234), QPointF(762.427, 10016188), QPointF(763.438, 4882192), QPointF(764.44, 1260325),
            QPointF(772.419, 1508312), QPointF(773.424, 420333), QPointF(1059.994, 267894), QPointF(1420.264, 327549)
    };

};

ScanPoints MS2ChargeDeconvolvotronTests::buildToyScanPoints() {

    const QVector<double> mzVals = BiophysicalCalcs::calculateIsotopesFromMz(1000.0, 2, 0, 3);
    const QVector<double> intzVals = {0.92, 1.0, 0.62, 0.29};

    ScanPoints scanPoints;
    Err e = ParallelUtils::zip(
            mzVals,
            intzVals,
            &scanPoints
    );

    return scanPoints;
}

void MS2ChargeDeconvolvotronTests::initTest() {

    ERR_INIT

    const QString &chargeModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

    const QString &monoModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

    const double ppmTol = 150.0;


    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
    e = ms2ChargeDeconvolvotron.init("KalliopeAndBellatrix4eva.dummyFile", "ChauncyAndFlops.dummpyfile",ppmTol);
    QCOMPARE(e, eFileError);

    e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);
    QCOMPARE(e, eNoError);

}

void MS2ChargeDeconvolvotronTests::buildMzValsForChargeMonoTest() {

    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;

    const QVector<double> mzValsExpexted = {
            998.997, 1000, 1001, 1002.01, 998.997, 999.498, 1000, 1000.5, 1001, 1001.51, 998.997, 999.331, 999.666, 1000, 1000.33, 1000.67, 1001, 1001.34
    };

    const QVector<double> mzVals = ms2ChargeDeconvolvotron.buildMzValsForChargeMono(1000.0);

    QCOMPARE(mzVals.size(), mzValsExpexted.size());

    for (int i = 0; i < mzVals.size(); i++) {
        QVERIFY(MathUtils::tSame(mzVals.at(i), mzValsExpexted.at(i)));
    }


}

void MS2ChargeDeconvolvotronTests::predictMS2ChargeAndMonoOffsetTest() {

    ERR_INIT

    const QString &chargeModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

    const QString &monoModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

    const double ppmTol = 15.0;


    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
    e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);
    QCOMPARE(e, eNoError);

    const ScanPoints scanPoints = buildToyScanPoints();

    const double mzExtract = 1000.5;

    int charge;
    MsUtils::chargeDeterminator({mzExtract, 1}, scanPoints, ppmTol, 1, 3, &charge);
    qDebug() << charge;

    const QVector<float> expectedPredVec = {4.09979e-07, 0.999999, 1.04794e-05, 0.0273309, 0.971944, 6.28135e-15};

    QVector<float> predVecCharge;
    QVector<float> predVecMono;
    e = ms2ChargeDeconvolvotron.testChargeMonoCaller(
            scanPoints,
            mzExtract,
            &predVecCharge,
            &predVecMono
            );
    QCOMPARE(e, eNoError);
    qDebug() << predVecCharge;
    qDebug() << predVecMono;

//    for (int i = 0; i < 3; i++) {
//        QVERIFY(MathUtils::tSame(predVec.at(i), expectedPredVec.at(i)));
//    }
}

void MS2ChargeDeconvolvotronTests::deisotopeScanPointsTest() {

    ERR_INIT

    const QString &chargeModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

    const QString &monoModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

    const double ppmTol = 15.0;

    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
    e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);
    QCOMPARE(e, eNoError);

    const ScanPoints scanPoints = buildToyScanPoints();

    ScanPoints scanPointsDeisotoped;
    e = ms2ChargeDeconvolvotron.deisotopeScanPoints(
            scanPoints,
            &scanPointsDeisotoped
            );
    QCOMPARE(e, eNoError);

    QVERIFY(MathUtils::tSame(scanPointsDeisotoped.front().x(), 1000.0));
    QVERIFY(MathUtils::tSame(scanPointsDeisotoped.front().y(), 0.9));

}

void MS2ChargeDeconvolvotronTests::deisotopeScanPointsRealDataTest() {

    ERR_INIT

    const QString &chargeModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

    const QString &monoModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

    const double ppmTol = 30.0;

    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
    e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);
    QCOMPARE(e, eNoError);

    ScanPoints scanPointsDeisotoped;
    e = ms2ChargeDeconvolvotron.deisotopeScanPoints(
            m_scanPoints,
            &scanPointsDeisotoped
    );
    QCOMPARE(e, eNoError);
//
//    e = MsUtils::writePointsToCSV(scanPointsDeisotoped, "/home/anichols/Desktop/deiso.csv");
//    QCOMPARE(e, eNoError);
//
//    e = MsUtils::writePointsToCSV(m_scanPoints, "/home/anichols/Desktop/points.csv");
//    QCOMPARE(e, eNoError);

    QCOMPARE(scanPointsDeisotoped.size(), 271);
    QCOMPARE(m_scanPoints.size(), 352);

}

QTEST_MAIN(MS2ChargeDeconvolvotronTests)

#include "MS2ChargeDeconvolvotronTests.moc"
