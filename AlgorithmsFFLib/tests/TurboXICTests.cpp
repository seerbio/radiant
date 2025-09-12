//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
// #include "MsFrame.h"
#include "MsReaderMzML.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"
#include "MsReaderMzMLLazyLoad.h"

#include <QtTest/QtTest>
#include <random>

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
	void readWriteTest();
	static void troubleShootTest();


private:

    QVector<MsScan> buildMsScans();

	QVector<MsScanInfo> m_msScanInfos;

};

QVector<MsScan> TurboXICTests::buildMsScans() {

	MsScanInfo msScanInfo;


    QVector<MsScan> msScans;

	MsScan msScan1;
	msScan1.mzVals = {100.11, 200.11};
	msScan1.intensityVals = {100.1, 200.1};
	msScanInfo.scanNumber = 1;
	msScanInfo.pointCount = 2;
	m_msScanInfos.push_back(msScanInfo);
	msScan1.msScanInfoPntr = &m_msScanInfos.back();
	msScans.push_back(msScan1);

	MsScan msScan2;
	msScan2.mzVals = {100.12, 200.12};
	msScan2.intensityVals = {100.2, 200.2};
	msScanInfo.scanNumber = 2;
	m_msScanInfos.push_back(msScanInfo);
	msScan2.msScanInfoPntr = &m_msScanInfos.back();
	msScans.push_back(msScan2);

	MsScan msScan3;
	msScan3.mzVals = {100.13, 200.13};
	msScan3.intensityVals = {100.3, 200.3};
	msScanInfo.scanNumber = 3;
	m_msScanInfos.push_back(msScanInfo);
	msScan3.msScanInfoPntr = &m_msScanInfos.back();
	msScans.push_back(msScan3);

	MsScan msScan4;
	msScan4.mzVals = {100.14, 200.14};
	msScan4.intensityVals = {100.4, 200.4};
	msScanInfo.scanNumber = 4;
	m_msScanInfos.push_back(msScanInfo);
	msScan4.msScanInfoPntr = &m_msScanInfos.back();
	msScans.push_back(msScan4);

	MsScan msScan5;
	msScan5.mzVals = {100.15, 200.15};
	msScan5.intensityVals = {100.5, 200.5};
	msScanInfo.scanNumber = 5;
	m_msScanInfos.push_back(msScanInfo);
	msScan5.msScanInfoPntr = &m_msScanInfos.back();
	msScans.push_back(msScan5);

    return msScans;
}

void TurboXICTests::initTest() {

    const QVector<MsScan> msScans = buildMsScans();

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(msScans);
    QCOMPARE(e, eNoError);

	const QVector<MsScan> msScansEmpty;
    e = turboXIC.init(msScansEmpty);
    QCOMPARE(e, eEmptyContainerError);

}

void TurboXICTests::extractPointsTest() {

    const QVector<MsScan> msScans = buildMsScans();

    ERR_INIT

    TurboXIC turboXIC;
    e = turboXIC.init(msScans);
    QCOMPARE(e, eNoError);

    XICPointsPntrs xicPoints;
	e = turboXIC.extractPointsXIC(100.0, 100.13, &xicPoints);
	QCOMPARE(e, eNoError);
    QCOMPARE(xicPoints.size(), 3);
    QCOMPARE(xicPoints.front()->scanNumber, 1);
    QVERIFY(MathUtils::tSame(xicPoints.front()->intensity, 100.1f));
    QVERIFY(MathUtils::tSame(xicPoints.back()->intensity, 100.3f));

}

void TurboXICTests::turboXICUtility() {

    ERR_INIT

    QSKIP("uncomment for troubleshooting");

//    const QString msDataFilePath
//            = QStringLiteral("/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq");
//
//    const MzTargetKey uniqueMsInfoScanKey = "454957";
//    double target = 454.957;
//    double window = 5.5;
//
//    MsFrame msFrame;
//    e = MsFrame::buildMsFrame(
//            msDataFilePath,
//            uniqueMsInfoScanKey,
//            {target - window, target + window},
//            &msFrame
//    );
//    QCOMPARE(e, eNoError);
//
//    e = msFrame.smoothFrame(
//            3,
//            1.0,
//            2,
//            1500.0
//            );
//    QCOMPARE(e, eNoError);
//
//    const QMap<int, QVector<PointFF>> &points = msFrame.scanNumberVsScanPoints();
//
//    TurboXIC turboXIC;
//    e = turboXIC.init(points);
//    QCOMPARE(e, eNoError);
//
//    const double mzCenter = 523.26232347;
//    const double ppmTol = 50;
//    const double massTol = MathUtils::calculatePPM(mzCenter, ppmTol);
//
//    const XICPoints xicPoints = turboXIC.extractPointsXIC(
//            mzCenter - massTol,
//            mzCenter + massTol,
//            0,
//            26000
//            );
//
//    qDebug() << xicPoints.scanNumbersVsIntensityVals.size();
//    const QVector<PointFF> vec = ParallelUtils::convertMapToPoints(xicPoints.scanNumbersVsIntensityVals);
//
//    for (const PointFF &p : vec) {
//        qDebug() << p;
//    }
//
//    e = MsUtils::writePointsToCSV(vec, "xic.csv");
//    QCOMPARE(e, eNoError);

}

void TurboXICTests::readWriteTest() {

	const QVector<MsScan> msScans = buildMsScans();

	ERR_INIT

	TurboXIC turboXIC;
	e = turboXIC.init(msScans);
	QCOMPARE(e, eNoError);

	XICPointsPntrs xicPoints;
	e = turboXIC.extractPointsXIC(100.0, 100.13, &xicPoints);
	QCOMPARE(e, eNoError);
	QCOMPARE(xicPoints.size(), 3);
	QCOMPARE(xicPoints.front()->scanNumber, 1);
	QVERIFY(MathUtils::tSame(xicPoints.front()->intensity, 100.1f));
	QVERIFY(MathUtils::tSame(xicPoints.back()->intensity, 100.3f));

	e = turboXIC.write("TurboXICReadWriteTest.xic");
	QCOMPARE(e, eNoError);

	TurboXIC turboXICRead;
	e = turboXICRead.init("TurboXICReadWriteTest.xic");
	QCOMPARE(e, eNoError);

	XICPointsPntrs xicPointsRead;
	e = turboXICRead.extractPointsXIC(100.0, 100.13, &xicPointsRead);
	QCOMPARE(e, eNoError);
	QCOMPARE(xicPointsRead.size(), 3);
	QCOMPARE(xicPointsRead.front()->scanNumber, 1);
	QVERIFY(MathUtils::tSame(xicPointsRead.front()->intensity, 100.1f));
	QVERIFY(MathUtils::tSame(xicPointsRead.back()->intensity, 100.3f));

}

namespace {

	Err openParallelTestLogic(
		const QString &fileName,
		const QVector<MsScanInfo*> &msScanInfos
		) {

		ERR_INIT

		QVector<MsScan> msScans;
		e = MsReaderMzMLLazyLoad::extractScanPoints(
			fileName,
			msScanInfos,
			&msScans
			); ree;

		TurboXIC turboXic;
		e = turboXic.init(msScans); ree;

		// std::random_device randomDevice;
		// std::mt19937 generator(randomDevice());
		// std::uniform_int_distribution<uint> distribution(0, std::numeric_limits<uint>::max());
		// const QString randomInt = QString::number(distribution(generator));
		//
		// const QString filePath = "TurboXICReadWriteTest_" + randomInt + ".xic";
		// e = turboXic.write(filePath); ree;
		//
		// e = turboXic.init(filePath); ree;

		qDebug() <<S_GLOBAL_TIMER.elapsed();

		ERR_RETURN
	}

	void deleteFilesWithExtension(const QString& directoryPath, const QString& extension) {
		QDir dir(directoryPath);

		if (!dir.exists()) {
			qDebug() << "Directory does not exist:" << directoryPath;
			return;
		}

		QStringList filters;
		filters << QString("*.%1").arg(extension); // e.g., *.xic
		QStringList files = dir.entryList(filters, QDir::Files);

		// Iterate through each file and delete it
		for (const QString& file : files) {
			QString filePath = dir.filePath(file); // Get the full path of the file
			QFile::remove(filePath);
			// ? qDebug() << "Deleted:" << filePath :
			//   qDebug() << "Failed to delete:" << filePath;
		}
	}

}//namespace
void TurboXICTests::troubleShootTest() {

	ERR_INIT

	QSKIP("Troubleshooting");

	const QString filename  = "/home/andrewnichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";

	MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
	e = msReaderMzMLLazyLoad.openFile(filename);
	QCOMPARE(e, eNoError);

	const QVector<MsScanInfo*> scanInfos = msReaderMzMLLazyLoad.getUniqueTandemMsScanInfos();

	QMap<MzTargetKey, QVector<MsScanInfo*>> mzTargetVsScanInfosPntrs = msReaderMzMLLazyLoad.m_mzTargetVsScanInfosPntrs;

	const auto binderLogic = std::bind(
		openParallelTestLogic,
		filename,
		std::placeholders::_1
		);

	QElapsedTimer et;
	et.start();
	QFuture<Err> futures = QtConcurrent::mapped(
		mzTargetVsScanInfosPntrs,
		binderLogic
		);
	futures.waitForFinished();
	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << et.elapsed() << "mSec to complete";

	deleteFilesWithExtension("/home/andrewnichols/Repos/Builds/PythiaDIA/bin", "xic");
}


QTEST_MAIN(TurboXICTests)
#include "TurboXICTests.moc"
