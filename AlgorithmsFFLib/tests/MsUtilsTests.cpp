//
// Created by anichols on 11/07/2021.
//

#include "MsUtils.h"
#include "ParallelUtils.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>

class MsUtilsTests : public QObject
{

Q_OBJECT

public:
    MsUtilsTests();
    ~MsUtilsTests() override = default;


private Q_SLOTS:

    void extractPointsFromPointsTest();
    void extractPointsFromPointsTestBST();
    static void extractPointsFromPointsSimpleTest();

private:

    const QVector<QPointF> m_points = {
            QPointF(101.071,71236.4), QPointF(101.108,17270.1), QPointF(101.253,2900.42), QPointF(102.055,63054.1),
            QPointF(103.055,3523.47), QPointF(103.059,4803.5), QPointF(104.053,5676.43), QPointF(106.276,40251.9),
            QPointF(106.387,12661.7), QPointF(107.357,2843.57), QPointF(110.072,47397.7), QPointF(111.075,4908.84),
            QPointF(112.087,3678), QPointF(115.087,12495.3), QPointF(115.123,3729.47), QPointF(116.071,5094), QPointF(120.081,547965),
            QPointF(121.078,3256.93), QPointF(121.085,46668), QPointF(122.072,3521.19), QPointF(123.511,3146.97),
            QPointF(127.051,3437.5), QPointF(127.087,4029.94), QPointF(129.066,5972.56), QPointF(129.103,434940),
            QPointF(129.105,5234.68), QPointF(130.05,8853.62), QPointF(130.084,6571.23), QPointF(130.087,289118),
            QPointF(130.106,16368.4), QPointF(131.09,10312.3), QPointF(131.118,5317.52), QPointF(133.043,8128.34),
            QPointF(136.076,32881.4), QPointF(139.087,4774.69), QPointF(141.066,5668.38), QPointF(143.118,139471),
            QPointF(144.122,10649.3), QPointF(145.061,5489.77), QPointF(147.113,579537), QPointF(148.076,9644.54),
            QPointF(148.111,4500.27), QPointF(148.117,31872.9), QPointF(152.082,7687.46), QPointF(155.082,24438.9),
            QPointF(157.098,7012.63), QPointF(157.134,5288.12), QPointF(158.093,9785.23), QPointF(159.077,20190.1),
            QPointF(159.092,11100.4), QPointF(165.103,9853.48), QPointF(166.062,11883.6), QPointF(166.087,3855.78),
            QPointF(167.093,12211.2), QPointF(169.098,9689.69), QPointF(170.061,4879.37), QPointF(170.167,3256.25),
            QPointF(171.077,9016.09), QPointF(171.113,32315.4), QPointF(171.15,3836.27), QPointF(172.043,8509.94),
            QPointF(173.093,77423.6), QPointF(173.129,10011.3), QPointF(174.096,10159.7), QPointF(174.132,3226.44),
            QPointF(175.087,3402.04), QPointF(175.12,46222.2), QPointF(177.103,7753.1), QPointF(179.093,8069.74),
            QPointF(180.077,5655.98), QPointF(181.061,5887.76), QPointF(181.109,4593.02), QPointF(181.131,2807.65),
            QPointF(183.077,3736.72), QPointF(183.114,4993.27), QPointF(183.15,4561.58), QPointF(185.093,4385.5),
            QPointF(185.129,6458.83), QPointF(186.092,3053.32), QPointF(186.124,5754.16), QPointF(187.072,40117.6),
            QPointF(187.108,11598.1), QPointF(188.075,3457.62), QPointF(188.14,3302), QPointF(191.119,46063.8),
            QPointF(192.122,3233.79), QPointF(193.098,7855.99), QPointF(195.088,12225), QPointF(197.104,28599.6),
            QPointF(197.129,12052.6), QPointF(197.166,12030.3), QPointF(198.089,3340.36), QPointF(199.054,5282.19),
            QPointF(199.072,13334.3), QPointF(199.181,5481.09), QPointF(200.103,4115.63), QPointF(201.088,34434.5),
            QPointF(201.124,24886.5), QPointF(202.091,9836.22), QPointF(205.097,5291.25), QPointF(207.088,10007),
            QPointF(208.072,23455.6), QPointF(209.057,3357.9), QPointF(209.104,13022.8), QPointF(212.104,5015.91),
            QPointF(213.088,6021.4), QPointF(214.12,4859.38), QPointF(215.103,9626.28), QPointF(215.14,4257.47),
            QPointF(216.107,2892.79), QPointF(216.135,4315.81), QPointF(217.065,26134.8), QPointF(217.083,23034.1),
            QPointF(217.098,16692), QPointF(217.12,10389.4), QPointF(217.134,1.40305e+06), QPointF(217.139,27167.4),
            QPointF(217.152,4906.18), QPointF(218.068,3118.85), QPointF(218.132,6574.8), QPointF(218.138,127723),
            QPointF(219.114,46562.6), QPointF(220.117,5861.07), QPointF(221.129,16978.5), QPointF(225.099,149961),
            QPointF(225.16,10354.2), QPointF(226.083,5356.45), QPointF(226.102,14783.9), QPointF(226.119,4401.36),
            QPointF(226.164,3151.3), QPointF(227.049,18529.6), QPointF(227.067,18384.2), QPointF(227.176,4936.14),
            QPointF(228.135,6303), QPointF(229.119,11710.8), QPointF(231.114,35572.2), QPointF(235.109,25473),
            QPointF(235.12,7701.21), QPointF(235.156,2871.48), QPointF(235.233,2882.39), QPointF(237.098,3162.32),
            QPointF(237.124,4198.71), QPointF(240.099,3245.16), QPointF(240.135,4901.91), QPointF(242.151,6800.69),
            QPointF(242.187,80888.8), QPointF(243.109,3417.56), QPointF(243.171,21927.6), QPointF(243.191,7240.04),
            QPointF(244.094,12614.4), QPointF(244.13,6471.31), QPointF(244.167,4561.52), QPointF(245.06,4095.64),
            QPointF(245.078,5181.47), QPointF(245.129,242110), QPointF(246.133,21263.8), QPointF(249.124,40507.9),
            QPointF(250.127,3651.83), QPointF(251.104,4248.86), QPointF(255.134,6789.11), QPointF(258.146,8873.95),
            QPointF(258.709,2857.99), QPointF(259.109,4731.34), QPointF(259.15,4771.3), QPointF(260.198,775939),
            QPointF(260.216,4503.08), QPointF(261.161,4980.26), QPointF(261.201,57549.8), QPointF(262.141,29143.5),
            QPointF(262.155,3588.65), QPointF(263.104,76231.9), QPointF(264.107,6763.61), QPointF(268.167,5340.88),
            QPointF(270.109,4449.07), QPointF(276.135,8282.35), QPointF(276.156,4035.77), QPointF(277.139,6817.35),
            QPointF(278.126,6600.08), QPointF(284.125,5139.9), QPointF(288.135,13688.8), QPointF(288.172,26500.2),
            QPointF(289.119,34496.8), QPointF(289.163,3109.68), QPointF(289.174,3431.7), QPointF(290.122,3413.23),
            QPointF(294.146,3168.17), QPointF(294.245,2666.09), QPointF(296.136,50694.9), QPointF(297.14,6549.51),
            QPointF(300.156,10974), QPointF(302.151,12711.7), QPointF(302.172,4378.64), QPointF(303.152,3384.44),
            QPointF(306.146,4531.16), QPointF(311.137,3142.61), QPointF(312.156,3876.21), QPointF(312.192,4324.11),
            QPointF(313.188,6270.01), QPointF(316.131,12945.6), QPointF(317.114,4953.58), QPointF(319.141,9707.28),
            QPointF(320.125,6891.53), QPointF(320.162,18716.9), QPointF(321.165,5330.5), QPointF(328.166,3152.61),
            QPointF(328.179,2764.98), QPointF(329.096,2897.45), QPointF(329.146,3063.56), QPointF(329.219,5399.85),
            QPointF(330.167,26712.3), QPointF(331.154,29541.7), QPointF(331.165,4722.37), QPointF(331.655,8461.5),
            QPointF(332.161,3730.14), QPointF(334.141,256925), QPointF(334.154,5250.42), QPointF(335.144,25702.7),
            QPointF(338.183,4762.7), QPointF(339.167,6306.63), QPointF(342.145,3567.34), QPointF(346.107,3014.96),
            QPointF(346.141,65101.8), QPointF(347.144,8524.23), QPointF(348.157,3731.47), QPointF(348.178,4236.41),
            QPointF(352.153,4241.17), QPointF(353.182,3550.96), QPointF(356.092,8216.62), QPointF(357.178,3090.66),
            QPointF(357.214,18366.2), QPointF(358.198,6791.22), QPointF(358.217,3674.21), QPointF(359.14,2923.41),
            QPointF(360.156,461319), QPointF(361.16,60021.5), QPointF(361.209,7661.25), QPointF(362.135,12537.4),
            QPointF(364.151,137975), QPointF(364.163,2270.73), QPointF(365.156,20416.7), QPointF(367.126,3373.27),
            QPointF(368.194,4071.47), QPointF(374.103,4101.51), QPointF(375.225,182115), QPointF(376.188,6092.25),
            QPointF(376.229,20721.5), QPointF(377.183,3791.72), QPointF(378.183,2737.91), QPointF(382.64,3002.97),
            QPointF(385.136,6358.34), QPointF(385.151,3897.21), QPointF(387.098,3993.45), QPointF(388.667,6102.52),
            QPointF(389.183,6928.14), QPointF(389.24,4484.56), QPointF(390.167,5507.83), QPointF(390.244,6376.11),
            QPointF(391.163,11116.8), QPointF(395.205,57540.9), QPointF(396.192,3228.37), QPointF(396.208,10575.7),
            QPointF(396.224,2836.67), QPointF(399.167,9496.28), QPointF(401.204,8497.99), QPointF(402.204,12506.9),
            QPointF(403.147,5110.78), QPointF(403.198,40892.8), QPointF(404.204,6494.47), QPointF(409.22,9533.7),
            QPointF(413.114,5341.25), QPointF(413.182,5884.1), QPointF(416.184,3611.15), QPointF(417.178,84420.6),
            QPointF(417.196,3624.68), QPointF(418.183,14176.9), QPointF(418.209,3360.16), QPointF(426.167,5093.79),
            QPointF(426.221,6333.34), QPointF(427.163,4790.51), QPointF(427.182,3331.49), QPointF(430.174,3206.91),
            QPointF(431.124,11398.9), QPointF(431.194,217524), QPointF(432.198,33451.9), QPointF(432.247,3952.26),
            QPointF(432.283,6776.78), QPointF(432.679,34299.9), QPointF(433.18,10548.9), QPointF(433.682,4623.9),
            QPointF(434.176,3372.66), QPointF(434.219,4437.76), QPointF(435.189,118600), QPointF(435.205,11554.5),
            QPointF(436.191,25229.7), QPointF(436.206,5434.4), QPointF(440.252,5195.63), QPointF(445.173,12103.7),
            QPointF(445.193,7128.48), QPointF(447.225,11333.2), QPointF(449.169,3702.62), QPointF(449.24,20904.9),
            QPointF(453.214,35880.8), QPointF(454.218,9050.84), QPointF(458.226,3689.81), QPointF(458.262,34424.1),
            QPointF(459.225,8413.91), QPointF(459.265,6685.69), QPointF(460.278,7698.06), QPointF(461.168,6248.12),
            QPointF(461.204,7671.58), QPointF(461.281,5377.93), QPointF(462.202,6998.98), QPointF(463.2,9162.33),
            QPointF(468.727,3738.22), QPointF(469.736,3846.75), QPointF(475.184,8158.42), QPointF(475.253,18360.5),
            QPointF(476.273,256857), QPointF(477.236,25427.9), QPointF(477.276,37016.5), QPointF(478.228,11866.6),
            QPointF(479.23,3655.9), QPointF(481.21,33319.2), QPointF(482.212,7153.65), QPointF(483.269,22722.7),
            QPointF(483.771,8708.46), QPointF(486.2,4820.26), QPointF(486.258,5237.15), QPointF(487.22,4654.39),
            QPointF(493.209,5300.71), QPointF(498.221,4281.37), QPointF(500.262,9791.43), QPointF(500.762,8999.12),
            QPointF(502.267,6125.6), QPointF(503.251,9610.09), QPointF(504.21,4937.16), QPointF(504.269,5206.42),
            QPointF(507.225,4070.83), QPointF(508.229,5005.97), QPointF(510.232,22558.2), QPointF(511.236,3404.58),
            QPointF(514.195,18974.3), QPointF(514.232,5654.47), QPointF(515.199,3687.98), QPointF(516.232,13457.7),
            QPointF(516.291,14958.4), QPointF(516.313,4189.87), QPointF(518.262,7781.11), QPointF(518.767,4073.52),
            QPointF(520.278,22066.4), QPointF(521.281,6791.3), QPointF(522.219,4022.19), QPointF(526.272,5116.43),
            QPointF(526.763,6031.89), QPointF(526.786,10407.9), QPointF(527.267,95189), QPointF(527.292,7217.1),
            QPointF(527.318,5114.9), QPointF(527.354,6829.98), QPointF(527.769,55433.2), QPointF(528.247,4145.07),
            QPointF(528.269,5240.66), QPointF(530.263,12121.3), QPointF(532.206,88859.7), QPointF(532.257,9369.76),
            QPointF(533.208,19557.4), QPointF(533.239,4871.46), QPointF(537.236,5462.47), QPointF(545.368,5611.69),
            QPointF(547.303,4602.96), QPointF(548.273,18234.2), QPointF(549.272,9770.44), QPointF(550.216,128242),
            QPointF(550.268,66324.8), QPointF(551.219,22675.9), QPointF(551.271,19712.1), QPointF(554.261,3851.12),
            QPointF(558.219,6696.3), QPointF(558.259,4171.02), QPointF(560.2,6457.54), QPointF(560.252,6018.68),
            QPointF(562.252,6558.48), QPointF(564.247,18296.5), QPointF(574.251,6579.21), QPointF(575.306,15827.7),
            QPointF(576.231,4281.67), QPointF(576.309,16819.9), QPointF(578.236,6124.56), QPointF(578.262,62657.9),
            QPointF(579.265,14260.6), QPointF(582.257,15540.3), QPointF(591.312,22053.8), QPointF(592.263,8165.96),
            QPointF(596.352,3539.43), QPointF(599.267,5449.28), QPointF(605.332,23932.5), QPointF(606.334,5252.68),
            QPointF(617.296,5539.06), QPointF(618.279,9101.79), QPointF(621.264,13367.9), QPointF(622.267,4631.29),
            QPointF(622.318,3902.83), QPointF(623.342,193157), QPointF(624.345,49996.2), QPointF(626.337,4715.83),
            QPointF(627.279,6864.19), QPointF(627.323,5658.06), QPointF(628.327,3165.98), QPointF(629.258,5291.35),
            QPointF(632.31,6823.52), QPointF(634.253,3087.34), QPointF(635.305,20668.7), QPointF(639.275,30630.1),
            QPointF(640.279,9130.23), QPointF(643.286,4032.95), QPointF(644.35,32220.1), QPointF(645.289,30005.3),
            QPointF(645.352,9877.48), QPointF(646.293,8486.67), QPointF(647.269,12444.3), QPointF(648.267,4311.67),
            QPointF(650.319,6756.76), QPointF(650.354,4869.9), QPointF(651.278,7587.69), QPointF(651.314,4036.45),
            QPointF(652.261,7144.93), QPointF(652.315,4028.9), QPointF(653.26,5040.95), QPointF(658.451,3866.67),
            QPointF(659.453,3680.16), QPointF(661.264,12654.6), QPointF(661.3,31320.5), QPointF(662.3,7211.28), QPointF(663.3,51532),
            QPointF(664.305,12237.4), QPointF(669.289,6636.31), QPointF(676.369,30665.9), QPointF(677.373,4440.68),
            QPointF(679.276,37310.4), QPointF(679.31,12558.8), QPointF(679.342,19950.7), QPointF(680.278,9031.21),
            QPointF(680.313,5910.28), QPointF(680.347,7733.38), QPointF(686.337,5872.56), QPointF(687.34,11695.9),
            QPointF(692.36,15734.3), QPointF(693.291,6990.31), QPointF(693.361,6573.51), QPointF(694.297,3567.02),
            QPointF(694.379,331685), QPointF(695.382,96347.8), QPointF(697.285,89019.5), QPointF(698.288,24976.4),
            QPointF(703.302,5612.35), QPointF(703.363,7124.83), QPointF(704.366,10176.2), QPointF(704.401,3552.49),
            QPointF(705.368,4438.92), QPointF(707.27,6071.52), QPointF(712.325,4950.14), QPointF(721.319,5155.96),
            QPointF(722.354,6907.17), QPointF(731.381,3074.4), QPointF(739.332,3604.34), QPointF(750.311,4877.37),
            QPointF(754.302,69204.6), QPointF(755.305,25890.4), QPointF(756.306,5221.11), QPointF(757.434,64615.5),
            QPointF(758.317,6103.96), QPointF(758.437,21104), QPointF(759.437,3469.68), QPointF(764.363,10327.9),
            QPointF(765.346,8720.07), QPointF(769.483,3224.56), QPointF(773.386,3471.42), QPointF(774.349,3728.67),
            QPointF(776.328,42746.1), QPointF(776.372,4855.73), QPointF(777.329,11350.7), QPointF(782.374,14799.4),
            QPointF(783.377,6571.16), QPointF(790.374,4601.53), QPointF(791.396,53259.7), QPointF(792.358,12443.8),
            QPointF(792.404,14172.7), QPointF(793.358,6333.01), QPointF(794.338,132471), QPointF(795.341,41622.6),
            QPointF(795.447,5008.72), QPointF(808.384,4714.58), QPointF(809.307,5249.35), QPointF(809.406,839768),
            QPointF(810.37,21854.7), QPointF(810.409,275910), QPointF(810.449,5296.91), QPointF(811.371,10426.4),
            QPointF(811.411,6211.28), QPointF(811.45,4415.47), QPointF(814.461,4940.38), QPointF(815.465,4120.32),
            QPointF(819.388,9296.57), QPointF(820.394,5282.53), QPointF(821.418,4060.89), QPointF(825.415,4653.38),
            QPointF(825.467,3524.7), QPointF(828.388,6655.34), QPointF(828.471,122453), QPointF(829.474,49792.7),
            QPointF(830.476,12100.1), QPointF(839.427,18329), QPointF(840.428,4726.12), QPointF(850.433,8769.49),
            QPointF(851.435,35792.6), QPointF(851.48,3526.87), QPointF(852.437,12822.3), QPointF(864.35,10549.2),
            QPointF(865.333,30157.8), QPointF(866.335,11976.7), QPointF(877.452,3281.46), QPointF(879.427,7420.95),
            QPointF(882.361,41699.4), QPointF(883.364,18912.5), QPointF(907.425,5638.05), QPointF(912.449,6987.76),
            QPointF(920.397,3007.22), QPointF(920.454,6397.92), QPointF(938.465,58945.3), QPointF(939.468,19539.2),
            QPointF(954.456,11677.1), QPointF(955.46,6147.91), QPointF(956.475,813249), QPointF(957.478,316979),
            QPointF(957.525,4151.42), QPointF(957.545,3964.87), QPointF(958.48,15241.1), QPointF(966.46,6803.66),
            QPointF(1008.17,3554.56), QPointF(1047.93,3332.15)
    };

    QVector<double> m_byIonsForPeptide_PFDAFTDLK = {
            360.155, 431.193, 578.261, 679.309, 794.336, 907.42, 1053.53, 375.224, 476.271, 623.34, 694.377, 809.404,
            956.472
    };

};

MsUtilsTests::MsUtilsTests() : QObject(){}

void MsUtilsTests::extractPointsFromPointsTest() {

    ERR_INIT

    const QVector<double> expectedYVals
            = {461319, 182115, 217524, 256857, 62657.9, 193157, 12558.8, 331685, 132471, 839768, 5638.05, 813249, -1};

    std::sort(m_byIonsForPeptide_PFDAFTDLK.begin(), m_byIonsForPeptide_PFDAFTDLK.end());

    QVector<QPointF> extractPoints;
    e = ParallelUtils::zip(m_byIonsForPeptide_PFDAFTDLK, expectedYVals, &extractPoints);
    QCOMPARE(e, eNoError);

	QElapsedTimer et;
	et.start();
    const ExtractPoints extractPointsVector
        = MsUtils::extractPointsFromPoints(m_points, extractPoints, 6.0);
	qDebug() << "Walking Extract elapsed time: " << et.nsecsElapsed() << "nSec";

    for (const QPointF &p : extractPointsVector.intensityFoundVsSearched) {
        QCOMPARE(p.x(), p.y());
    }

}

void MsUtilsTests::extractPointsFromPointsTestBST() {

	ERR_INIT

	const QVector<double> expectedYVals
			= {461319, 182115, 217524, 256857, 62657.9, 193157, 12558.8, 331685, 132471, 839768, 5638.05, 813249, -1};

	std::sort(m_byIonsForPeptide_PFDAFTDLK.begin(), m_byIonsForPeptide_PFDAFTDLK.end());

	QVector<QPointF> extractPoints;
	e = ParallelUtils::zip(m_byIonsForPeptide_PFDAFTDLK, expectedYVals, &extractPoints);
	QCOMPARE(e, eNoError);

	QElapsedTimer et;
	et.start();
	const ExtractPoints extractPointsVector
		= MsUtils::extractPointsFromPointsBST(m_points, extractPoints, 6.0);
	qDebug() << "BST Extract elapsed time: " << et.nsecsElapsed() << "nSec";

	for (const QPointF &p : extractPointsVector.intensityFoundVsSearched) {
		QCOMPARE(p.x(), p.y());
	}

}

void MsUtilsTests::extractPointsFromPointsSimpleTest() {

    QVector<QPointF> p1 = {{1, 10}, {2,20}, {3,30}};
    QVector<QPointF> p2 = {{1, 10}, {2,20}, {4,40}};

    const ExtractPoints extractPointsVector
            = MsUtils::extractPointsFromPoints(p1, p2, 10.0);

    const QVector<QPointF> expectedXVals = {QPointF(1,1), QPointF(2,2), QPointF(-1,4)};
    const QVector<QPointF> expectedYVals = {QPointF(10,10), QPointF(20,20), QPointF(-1,40)};
    QCOMPARE(extractPointsVector.mzFoundVsSearched.size(), expectedXVals.size());
    QCOMPARE(extractPointsVector.intensityFoundVsSearched.size(), expectedYVals.size());

    for (int i = 0; i < extractPointsVector.mzFoundVsSearched.size(); i++) {
        QCOMPARE(extractPointsVector.mzFoundVsSearched.at(i), expectedXVals.at(i));
        QCOMPARE(extractPointsVector.intensityFoundVsSearched.at(i), expectedYVals.at(i));
    }
}


QTEST_MAIN(MsUtilsTests)
#include "MsUtilsTests.moc"
