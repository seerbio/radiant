//
// Created by anichols on 11/07/2021.
//

#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"

#include <QtTest/QtTest>


class MsReaderBaseTests : public QObject
{
    Q_OBJECT

public:
    MsReaderBaseTests() = default;
    ~MsReaderBaseTests() override = default;

private Q_SLOTS:

    static void setGetMsScanInfoTest();
    static void setGetScanPointsTest();
    static void resetTest();
    static void closeFileTest();
    void filePathTest();
    void isDIATest();
    void isInitTest();
    void scanTimeMinMaxTest();
    void getScanPointsTest();
    void getScanPointsPntrsTest();
    void getScanPointsMsLevelTest1();
    void getScanPointsMsLevelTest2();
    void getMsScanInfoTest();
    void collateMS2MzTargetFramesTest();
    void getUniqueTandemMsScanInfos();
    void getFrameCountTest();
    static void sortScanPointsTest();
    void getNearestScanNumberFromScanTimeTest();
    void getScanNumberVsScanTimeTest();
    static void splitScanPointsTest();
    static void zipScanPointsTest();


private:

    const QString m_prqFFFilePath
        = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

};

void MsReaderBaseTests::setGetMsScanInfoTest() {

    QMap<ScanNumber, MsScanInfo> msScanInfos = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    MsReaderBase msReaderBase;
    msReaderBase.setMsScanInfo(msScanInfos);

    QCOMPARE(msReaderBase.getMsScanInfos().size(), msScanInfos.size());

}

void MsReaderBaseTests::setGetScanPointsTest() {

    ERR_INIT

    QMap<ScanNumber, MsScanInfo> msScanInfos = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    MsReaderBase msReaderBase;
    msReaderBase.setMsScanInfo(msScanInfos);

    QCOMPARE(msReaderBase.getMsScanInfos().size(), msScanInfos.size());

    QMap<ScanNumber, ScanPoints> scanPoints = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    e = msReaderBase.setScanPoints(scanPoints);
    QCOMPARE(e, eNoError);

    QCOMPARE(msReaderBase.getScanPoints().size(), scanPoints.size());
}

void MsReaderBaseTests::resetTest() {

    QMap<ScanNumber, MsScanInfo> msScanInfos = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    QMap<ScanNumber, ScanPoints> scanPoints = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    MsReaderBase msReaderBase;
    msReaderBase.setMsScanInfo(msScanInfos);
    msReaderBase.setScanPoints(scanPoints);

    QCOMPARE(msReaderBase.getScanPoints().size(), scanPoints.size());
    QCOMPARE(msReaderBase.getMsScanInfos().size(), msScanInfos.size());

    msReaderBase.reset();
    QCOMPARE(msReaderBase.getScanPoints().size(), 0);
    QCOMPARE(msReaderBase.getMsScanInfos().size(), 0);

}

void MsReaderBaseTests::closeFileTest() {

    QMap<ScanNumber, MsScanInfo> msScanInfos = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    QMap<ScanNumber, ScanPoints> scanPoints = {
            {1, {}},
            {2, {}},
            {3, {}},
    };

    MsReaderBase msReaderBase;
    msReaderBase.setMsScanInfo(msScanInfos);
    msReaderBase.setScanPoints(scanPoints);

    QCOMPARE(msReaderBase.getScanPoints().size(), scanPoints.size());
    QCOMPARE(msReaderBase.getMsScanInfos().size(), msScanInfos.size());

    msReaderBase.closeFile();
    QCOMPARE(msReaderBase.getScanPoints().size(), 0);
    QCOMPARE(msReaderBase.getMsScanInfos().size(), 0);

}

void MsReaderBaseTests::filePathTest() {

    MsReaderPointerAcc msReaderPointerAcc;
    msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), m_prqFFFilePath);

}

void MsReaderBaseTests::isDIATest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const bool msFileIsDIA = msReaderPointerAcc.ptr->isDIA();
    QCOMPARE(msFileIsDIA, true);
}

void MsReaderBaseTests::isInitTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);
    QCOMPARE(msReaderPointerAcc.isInit(), true);

    MsReaderPointerAcc msReaderPointerAccFalse;
    QCOMPARE(msReaderPointerAccFalse.isInit(), false);

}

void MsReaderBaseTests::scanTimeMinMaxTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const QPair<ScanTime, ScanTime> scanTimeMinMax = msReaderPointerAcc.ptr->scanTimeMinMax();

    QCOMPARE(MathUtils::pRound(static_cast<double>(scanTimeMinMax.first), 4), 0.0028);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanTimeMinMax.second), 4), 6.1676);

}

void MsReaderBaseTests::getScanPointsPntrsTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs = msReaderPointerAcc.ptr->getScanPointsPntrs();
    QCOMPARE(scanNumberVsScanPointsPntrs.size(), 5000);
    QCOMPARE(scanNumberVsScanPointsPntrs.first()->size(), 49);
    QCOMPARE(scanNumberVsScanPointsPntrs.firstKey(), 1);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanPointsPntrs.first()->at(0).x()), 2), 380.95);
    QCOMPARE(static_cast<int>(scanNumberVsScanPointsPntrs.first()->at(0).y()), 124992);
}

void MsReaderBaseTests::getScanPointsTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints = msReaderPointerAcc.ptr->getScanPoints();
    QCOMPARE(scanNumberVsScanPoints.size(), 5000);
    QCOMPARE(scanNumberVsScanPoints.first().size(), 49);
    QCOMPARE(scanNumberVsScanPoints.firstKey(), 1);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanPoints.first().at(0).x()), 2), 380.95);
    QCOMPARE(static_cast<int>(scanNumberVsScanPoints.first().at(0).y()), 124992);
}

void MsReaderBaseTests::getScanPointsMsLevelTest1() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
    e = msReaderPointerAcc.ptr->getScanPoints(1, &scanNumberVsScanPoints);
    QCOMPARE(e, eNoError);
    QCOMPARE(scanNumberVsScanPoints.size(), 80);
    QCOMPARE(scanNumberVsScanPoints.first().size(), 49);
    QCOMPARE(scanNumberVsScanPoints.firstKey(), 1);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanPoints.first().at(0).x()), 2), 380.95);
    QCOMPARE(static_cast<int>(scanNumberVsScanPoints.first().at(0).y()), 124992);

}

void MsReaderBaseTests::getScanPointsMsLevelTest2() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs;
    e = msReaderPointerAcc.ptr->getScanPoints(1, &scanNumberVsScanPointsPntrs);
    QCOMPARE(e, eNoError);
    QCOMPARE(scanNumberVsScanPointsPntrs.size(), 80);
    QCOMPARE(scanNumberVsScanPointsPntrs.first()->size(), 49);
    QCOMPARE(scanNumberVsScanPointsPntrs.firstKey(), 1);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanPointsPntrs.first()->at(0).x()), 2), 380.95);
    QCOMPARE(static_cast<int>(scanNumberVsScanPointsPntrs.first()->at(0).y()), 124992);

}

void MsReaderBaseTests::getMsScanInfoTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfo;
    e = msReaderPointerAcc.ptr->getMsScanInfo(666, &msScanInfo);
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfo.scanNumber, 666);
    QCOMPARE(msScanInfo.msLevel, 2);
    QCOMPARE(MathUtils::pRound(static_cast<double>(msScanInfo.scanTime), 4), 0.8267);
    QCOMPARE(msScanInfo.collisionEnergy, 28);
    QCOMPARE(MathUtils::pRound(static_cast<double>(msScanInfo.precursorTargetMz), 3), 725.079);
    QCOMPARE(msScanInfo.isoWindowUpper, 5.5);
    QCOMPARE(msScanInfo.isoWindowLower, 5.5);
    QCOMPARE(msScanInfo.ionMobilityDriftTime, -1.0);
    QCOMPARE(msScanInfo.ionMobilityIndex, -1);
    QCOMPARE(msScanInfo.targetKey(), "725079");

}

void MsReaderBaseTests::collateMS2MzTargetFramesTest() {

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(&diaTargetFrames);
    QCOMPARE(e, eNoError);
    QCOMPARE(diaTargetFrames.size(), 62);

}

void MsReaderBaseTests::getUniqueTandemMsScanInfos() {

    ERR_INIT;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    QVector<MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getUniqueTandemMsScanInfos();
    QCOMPARE(msScanInfos.size(), 62);

}

void MsReaderBaseTests::getFrameCountTest() {

    ERR_INIT;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    QCOMPARE(msReaderPointerAcc.ptr->getFrameCount(), 62);
}

void MsReaderBaseTests::sortScanPointsTest() {

    ScanPoints scanPoints = {
            {100.0, 666.0},
            {200.0, 266.0},
            {300.0, 1266.0},
            {400.0, 66.0}
    };

    MsReaderBase::sortScanPoints(
            ScanPointsSort::AscIntensity,
            &scanPoints
            );
    QCOMPARE(scanPoints.front().y(), 66.0);
    QCOMPARE(scanPoints.back().y(), 1266.0);

    MsReaderBase::sortScanPoints(
            ScanPointsSort::DescIntensity,
            &scanPoints
    );
    QCOMPARE(scanPoints.back().y(), 66.0);
    QCOMPARE(scanPoints.front().y(), 1266.0);


    MsReaderBase::sortScanPoints(
            ScanPointsSort::AscMz,
            &scanPoints
    );
    QCOMPARE(scanPoints.front().x(), 100.0);
    QCOMPARE(scanPoints.back().x(), 400.0);

    MsReaderBase::sortScanPoints(
            ScanPointsSort::DescMz,
            &scanPoints
    );
    QCOMPARE(scanPoints.back().x(), 100.0);
    QCOMPARE(scanPoints.front().x(), 400.0);
}

void MsReaderBaseTests::getNearestScanNumberFromScanTimeTest() {

    ERR_INIT;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const ScanTime scanTime = 0.666;
    const ScanNumber scanNumber = msReaderPointerAcc.ptr->getNearestScanNumberFromScanTime(scanTime);
    QCOMPARE(scanNumber, 536);

    MsScanInfo msScanInfoBefore;
    e = msReaderPointerAcc.ptr->getMsScanInfo(scanNumber - 1, &msScanInfoBefore);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfoTarget;
    e = msReaderPointerAcc.ptr->getMsScanInfo(scanNumber, &msScanInfoTarget);
    QCOMPARE(e, eNoError);

    MsScanInfo msScanInfoAfter;
    e = msReaderPointerAcc.ptr->getMsScanInfo(scanNumber + 1, &msScanInfoAfter);
    QCOMPARE(e, eNoError);

    QVERIFY(std::abs(msScanInfoBefore.scanTime - scanTime) > std::abs(msScanInfoTarget.scanTime - scanTime));
    QVERIFY(std::abs(msScanInfoAfter.scanTime - scanTime) > std::abs(msScanInfoTarget.scanTime - scanTime));
}

void MsReaderBaseTests::getScanNumberVsScanTimeTest() {

    ERR_INIT;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFFFilePath);
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanTime> scanNumberVsScanTime = msReaderPointerAcc.ptr->getScanNumberVsScanTime();
    QCOMPARE(scanNumberVsScanTime.size(), 5000);
    QCOMPARE(scanNumberVsScanTime.firstKey(), 1);
    QCOMPARE(scanNumberVsScanTime.lastKey(), 5000);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanTime.first()), 4), 0.0028);
    QCOMPARE(MathUtils::pRound(static_cast<double>(scanNumberVsScanTime.last()), 4), 6.1676);

}

void MsReaderBaseTests::splitScanPointsTest() {

    ERR_INIT

    ScanPoints scanPoints = {
            {100.0, 666.0},
            {200.0, 266.0},
            {300.0, 1266.0},
            {400.0, 66.0}
    };

    QVector<float> mzVals;
    QVector<float> intensityVals;
    e = MsReaderBase::splitScanPoints(scanPoints, &mzVals, &intensityVals);
    QCOMPARE(e, eNoError);
    QCOMPARE(mzVals.size(), 4);
    QCOMPARE(intensityVals.size(), 4);
    QCOMPARE(scanPoints.front().x(), 100.0);
    QCOMPARE(scanPoints.back().x(), 400.0);

    e = MsReaderBase::splitScanPoints(&scanPoints, &mzVals, &intensityVals);
    QCOMPARE(e, eNoError);
    QCOMPARE(mzVals.size(), 4);
    QCOMPARE(intensityVals.size(), 4);
    QCOMPARE(scanPoints.front().x(), 100.0);
    QCOMPARE(scanPoints.back().x(), 400.0);

    scanPoints.clear();
    e = MsReaderBase::splitScanPoints(scanPoints, &mzVals, &intensityVals);
    QCOMPARE(e, eEmptyContainerError);

}

void MsReaderBaseTests::zipScanPointsTest() {

    ERR_INIT

    QVector<float> mzVals = {100.0, 200.0, 300.0};
    QVector<float> intnesityVals = {100.0, 200.0, 300.0};

    ScanPoints scanPoints;
    e = MsReaderBase::zipScanPoints(mzVals, intnesityVals, &scanPoints);
    QCOMPARE(e, eNoError);
    QCOMPARE(scanPoints.size(), mzVals.size());

    QCOMPARE(scanPoints.front().x(), mzVals.front());
    QCOMPARE(scanPoints.front().y(), intnesityVals.front());
    QCOMPARE(scanPoints.back().x(), mzVals.back());
    QCOMPARE(scanPoints.back().y(), intnesityVals.back());

    mzVals.pop_back();
    e = MsReaderBase::zipScanPoints(mzVals, intnesityVals, &scanPoints);
    QCOMPARE(e, eError);

}


QTEST_MAIN(MsReaderBaseTests)
#include "MsReaderBaseTests.moc"







