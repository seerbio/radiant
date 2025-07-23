//
// Created by anichols on 11/07/2021.
//

#include "MsReaderMzMLLazyLoad.h"
#include "FastaReader.h"
#include "ParallelUtils.h"

#include <QString>
#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>


class MsReaderMZMLLazyLoadTests : public QObject
{
    Q_OBJECT

public:
    MsReaderMZMLLazyLoadTests() = default;
    ~MsReaderMZMLLazyLoadTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    static void openFileTest();
    static void troubleShoot();

};


void MsReaderMZMLLazyLoadTests::openFileTest() {

    QSKIP("Turn on after complete and making tests.");

    // const QString &msParquetFilePath = QDir(qApp->applicationDirPath()).filePath("1min.mzML");
    //
    // ERR_INIT
    //
    // QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));
    //
    // MsReaderMzML reader;
    // e = reader.openFile(msParquetFilePath);
    // QCOMPARE(e, Error::eNoError);
    // QCOMPARE(reader.m_msScanInfo.size(), 372);
    // QCOMPARE(reader.m_scanPoints.size(), 372);
    //
    // MsScanInfo msScanInfo;
    // e = reader.getMsScanInfo(372, &msScanInfo);
    // QCOMPARE(e, eNoError);
    // QCOMPARE(msScanInfo.scanNumber, 372);
    // QCOMPARE(msScanInfo.msLevel, 1);
    // QCOMPARE(msScanInfo.collisionEnergy, -1);
    // QCOMPARE(QString::number(msScanInfo.scanTime), "10.0411");
    // QCOMPARE(QString::number(msScanInfo.precursorTargetMz), "-1");
    // QCOMPARE(msScanInfo.isoWindowLower, -1);
    // QCOMPARE(msScanInfo.isoWindowUpper, -1);
    //
    // QPair<Err, ScanPoints*> scanPointsResult = reader.getScanPoints(msScanInfo.scanNumber);
    // e = scanPointsResult.first;
    // QCOMPARE(e, eNoError);
    //
    // const ScanPoints scanPoints = *scanPointsResult.second;
    //
    // QCOMPARE(int(scanPoints.first().x()), 150);
    // QCOMPARE(int(scanPoints.first().y()), 814);
    // QCOMPARE(int(scanPoints.last().x()), 1250);
    // QCOMPARE(int(scanPoints.last().y()), 2646);
    // QCOMPARE(scanPoints.size(), 2092);

}

namespace {

	int openParallelTestLogic(
		const QString &fileName,
		const QVector<MsScanInfo*> &msScanInfos
		) {

    	ERR_INIT

    	QVector<MsScan> msScans;
    	e = MsReaderMzMLLazyLoad::extractScanPoints(
			fileName,
			msScanInfos,
			&msScans
			);

    	// qDebug() << e << msScans.size() << "scanNumberVsScanPoints Size yo!" << msScans.front().mzVals.at(0);
    	// qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "ms";

    	return 1;
    }

}//namespace
void MsReaderMZMLLazyLoadTests::troubleShoot() {

    ERR_INIT

    // QSKIP("Troubleshooting");
    //TODO make more tests.

    const QString filename  = "/home/ubuntu/Data/EXP23111_2023ms0979bX45_A.raw.mzML";
    // const QString filename = "/home/andrewnichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";
    // const QString filename = "/home/andrewnichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML";

    // const QString &filename = QDir(qApp->applicationDirPath()).filePath("1min.mzML");

    MsReaderMzMLLazyLoad msReaderMzMLLazyLoad;
    e = msReaderMzMLLazyLoad.openFile(filename);
    QCOMPARE(e, eNoError);
    // QCOMPARE(msReaderMzMLLazyLoad.m_mzTargetVsScanInfosPntrs.size(), 201);

	QMap<MzTargetKey, QVector<MsScanInfo*>> mzTargetVsScanInfosPntrs = msReaderMzMLLazyLoad.m_mzTargetVsScanInfosPntrs;

	const auto binderLogic = std::bind(
		openParallelTestLogic,
		filename,
		std::placeholders::_1
		);

	QElapsedTimer et;
	et.start();
	QFuture<int> futures = QtConcurrent::mapped(
		mzTargetVsScanInfosPntrs,
		binderLogic
		);
	futures.waitForFinished();
	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << et.elapsed();


}


QTEST_MAIN(MsReaderMZMLLazyLoadTests)
#include "MsReaderMZMLLazyLoadTests.moc"
