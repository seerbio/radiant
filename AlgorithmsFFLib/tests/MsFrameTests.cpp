#include "MsFrame.h"

#include "MsReaderParquet.h"

#include <QElapsedTimer>
#include <QVector>
#include <QtTest>


class MsFrameTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameTests() = default;
    ~MsFrameTests() override = default;


private slots:

    void initTest();
    void isValidTest();
    void writeFrameScansTest();
    void scanCountTest();
    void frameIndexVsScanPointsTest();
    void scanNumberVsScanPointsTest();
    void scanNumberFromFrameIndexTest();
    void scanTimeFromScanNumberTest();
    void scanNumberFromScanTimeTest();
    void frameIndexFromScanNumberTest();
    void getScanPointsByScanNumberTest();

private:

    QMap<ScanNumber, ScanPoints> m_scanPoints = {
            {1, {{100.0, 1000.0}, {200.0, 2000.0}, {300.0, 3000.0}}},
            {10, {{100.0, 1000.0}, {200.0, 2000.0}, {300.0, 3000.0}}},
            {20, {{101.0, 1001.0}, {201.0, 2001.0}, {301.0, 3001.0}}},
            {30, {{100.0, 1000.0}, {200.0, 2000.0}, {300.0, 3000.0}}}
    };

    const QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime = {
            {1, 10.0},
            {10, 20.0},
            {20, 30.0},
            {30, 40.0},
    };

};


void MsFrameTests::initTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;

    e = msFrame.init(
            scanPointsPtrs,
            {}
    );
    QCOMPARE(e, eEmptyContainerError);

    e = msFrame.init(
            {},
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eEmptyContainerError);

    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

}

void MsFrameTests::isValidTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;

    const bool isValidFail = msFrame.isValid();
    QCOMPARE(isValidFail, false);

    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

    const bool isValidPass = msFrame.isValid();
    QCOMPARE(isValidPass, true);

}

void MsFrameTests::writeFrameScansTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    e = MsFrame::writeFrameScans(m_scanPoints, "test1.frame");
    QCOMPARE(e, eNoError);

    e = ErrorUtils::fileExists("test1.frame");
    QCOMPARE(e, eNoError);

    e = MsFrame::writeFrameScans(scanPointsPtrs, "test2.frame");
    QCOMPARE(e, eNoError);

    e = ErrorUtils::fileExists("test2.frame");
    QCOMPARE(e, eNoError);

}

void MsFrameTests::scanCountTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(msFrame.scanCount(), scanPointsPtrs.size());

}

void MsFrameTests::frameIndexVsScanPointsTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

    const QMap<FrameIndex, ScanPoints*> frameIndexVsScanPoints = msFrame.frameIndexVsScanPoints();
    const QList<FrameIndex> expectedFrameIndexes = {0, 1, 2, 3};
    QCOMPARE(frameIndexVsScanPoints.keys(), expectedFrameIndexes);

}

void MsFrameTests::scanNumberVsScanPointsTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints = msFrame.scanNumberVsScanPoints();
    const QList<ScanNumber> expectedScanNumbers = {1, 10, 20, 30};
    QCOMPARE(scanNumberVsScanPoints.keys(), expectedScanNumbers);

}

void MsFrameTests::scanNumberFromFrameIndexTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(msFrame.scanNumberFromFrameIndex(0), 1);
    QCOMPARE(msFrame.scanNumberFromFrameIndex(1), 10);
    QCOMPARE(msFrame.scanNumberFromFrameIndex(2), 20);
    QCOMPARE(msFrame.scanNumberFromFrameIndex(3), 30);
}

void MsFrameTests::scanTimeFromScanNumberTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(msFrame.scanTimeFromScanNumber(20), 30);
    QCOMPARE(msFrame.scanTimeFromScanNumber(21), 0);
}

void MsFrameTests::scanNumberFromScanTimeTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(msFrame.scanNumberFromScanTime(30.1), 20);

}

void MsFrameTests::frameIndexFromScanNumberTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(msFrame.frameIndexFromScanNumber(20), 2);

}

void MsFrameTests::getScanPointsByScanNumberTest() {

    ERR_INIT

    QMap<ScanNumber, ScanPoints*> scanPointsPtrs;
    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrame;
    e = msFrame.init(
            scanPointsPtrs,
            m_scanNumberVsScanTime
    );
    QCOMPARE(e, eNoError);

    const ScanPoints* scanPointsPtrsReturn = msFrame.getScanPointsByScanNumber(20);
    QCOMPARE(scanPointsPtrsReturn->at(0).x(), 101.0);
    QCOMPARE(scanPointsPtrsReturn->at(0).y(), 1001.0);

}


QTEST_MAIN(MsFrameTests)

#include "MsFrameTests.moc"
