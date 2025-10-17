//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtTest/QtTest>

#include "MS2Ion.h"
#include "ObjectCSVWriters.h"

class TurboXICTests : public QObject
{
    Q_OBJECT

public:
    TurboXICTests() = default;
    ~TurboXICTests() override = default;

private Q_SLOTS:

    void initTest();
    void extractPointsTest();
    void turboXICUtility();


private:

    QMap<int, QVector<PointFF>> buildPoints() const;

};


QMap<int, QVector<PointFF>> TurboXICTests::buildPoints() const {

    QMap<int, QVector<PointFF>> points;
    points.insert(1, {PointFF(100.11, 100.1), PointFF(200.11, 200.1)});
    points.insert(2, {PointFF(100.12, 100.2), PointFF(200.12, 200.2)});
    points.insert(3, {PointFF(100.13, 100.3), PointFF(200.13, 200.3)});
    points.insert(4, {PointFF(100.14, 100.4), PointFF(200.14, 200.4)});
    points.insert(5, {PointFF(100.15, 100.5), PointFF(100.151, 100.5), PointFF(200.15, 200.5)});


    return points;
}

void TurboXICTests::initTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanPoints*> emptyPoints;
    e = turboXIC.init(emptyPoints);
    QCOMPARE(e, eEmptyContainerError);

}

void TurboXICTests::extractPointsTest() {

    QMap<int, QVector<PointFF>> _points = buildPoints();

    QMap<int, QVector<PointFF>*> points;
    for (auto it = _points.begin(); it != _points.end(); it++) {
        points.insert(it.key(), &it.value());
    }

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(points);
    QCOMPARE(e, eNoError);

    //TODO fix test

    const XICPoints xicPoints = turboXIC.extractPointsXIC(100.0, 100.13);

    QCOMPARE(xicPoints.size(), 3);
    QCOMPARE(xicPoints.front().scanNumber, 1);
    QVERIFY(MathUtils::tSame(xicPoints.front().intensity, 100.1f));
    QVERIFY(MathUtils::tSame(xicPoints.back().intensity, 100.3f));

}

void TurboXICTests::turboXICUtility() {

    ERR_INIT

	const QString filePath = "/home/andrewnichols/Desktop/Data/MsData/EXP24078_2024EU0002rX1_A.raw.mzML";
	MsReaderPointerAcc msReaderPointerAcc;
	e = msReaderPointerAcc.openFile(filePath);
	QCOMPARE(e, eNoError);

	QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
	e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(&diaTargetFrames);
	QCOMPARE(e, eNoError);

	// "480969" "KLVPFATELHER" false 3 -1 0 0 0 0 0 0 0 0 0 0 0 0 SDLFJDSL

	MzTargetKey targetKey = "480969";
	const QMap<ScanNumber, ScanPoints*> &scanPoints = diaTargetFrames.value(targetKey);


	MsFrame msFrame;
	e = msFrame.init(scanPoints, msReaderPointerAcc.ptr->getScanNumberVsScanTime());
	QCOMPARE(e, eNoError);

	TurboXIC turboXIC;
	e = turboXIC.init(msFrame.frameIndexVsScanPoints());
	QCOMPARE(e, eNoError);

	const QVector<float> mzVals = {
		550.28,
		599.814,
		441.22,
		341.255,
		554.305,
		784.395,
		855.432,
		683.347,
		277.656,
		392.701,
		656.356,
		1099.55
	};

	for (int i = 0; i < mzVals.size(); i++) {

		const float mzVal = mzVals.at(i);

		const float ppmTol = 11;
		const float mzTol = MathUtils::calculatePPM(mzVal, ppmTol);
		const XICPoints xic = turboXIC.extractPointsXIC(
			mzVal - mzTol,
			mzVal + mzTol
			);

		QVector<float> vec(1000, 0.0);
		for (const XICPoint &pt : xic) {
			vec[pt.scanNumber] += pt.intensity;
		}

		const QString outputPath = QStringLiteral("vec%1").arg(i);
		e = ObjectCSVWriters::writeVectorToFile(vec, outputPath);
		QCOMPARE(e, eNoError);
	}

	qDebug()
	<< msFrame.scanTimeFromFrameIndex(425)
	<< msFrame.scanTimeFromFrameIndex(433)
	;
}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
