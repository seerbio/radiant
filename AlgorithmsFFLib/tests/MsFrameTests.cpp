#include "MsFrame.h"
#include "MsFrameLookupReference.h"

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
    void cachedFrameLookupsMatchLegacy_data();
    void cachedFrameLookupsMatchLegacy();

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


void MsFrameTests::cachedFrameLookupsMatchLegacy_data() {
    QTest::addColumn<int>("count");
    QTest::addColumn<int>("kind");
    for (int count : {1, 2, 3, 15, 16, 31, 257, 1423}) {
        for (int kind = 0; kind < 9; ++kind) {
            const QByteArray name = QByteArray::number(count) + "-" + QByteArray::number(kind);
            QTest::newRow(name.constData()) << count << kind;
        }
    }
}

void MsFrameTests::cachedFrameLookupsMatchLegacy() {
    QFETCH(int, count);
    QFETCH(int, kind);
    MsFrame actual;
    MsFrameLookupReference expected;
    for (int initialization = 0; initialization < 2; ++initialization) {
        auto inputs = makeFrameLookupInputs(initialization ? count / 2 + 1 : count,
                                           initialization ? (kind + 4) % 9 : kind);
        QMap<ScanNumber, ScanPoints *> pointers;
        for (auto it = inputs.points.begin(); it != inputs.points.end(); ++it)
            pointers.insert(it.key(), &it.value());
        expected.init(inputs.points.keys(), inputs.times);
        QCOMPARE(actual.init(pointers, inputs.times), eNoError);
        QCOMPARE(actual.scanCount(), expected.scanCount());
        for (int frame : {std::numeric_limits<int>::min(), -10, -1, 0, 1,
                          expected.scanCount() - 1, expected.scanCount(),
                          expected.scanCount() + 1, std::numeric_limits<int>::max()}) {
            QCOMPARE(actual.scanNumberFromFrameIndex(frame), expected.scanNumberFromFrameIndex(frame));
            QCOMPARE(scanTimeBits(actual.scanTimeFromFrameIndex(frame)),
                     scanTimeBits(expected.scanTimeFromFrameIndex(frame)));
        }
        for (int frame = 0; frame < expected.scanCount(); ++frame) {
            QCOMPARE(actual.scanNumberFromFrameIndex(frame), expected.scanNumberFromFrameIndex(frame));
            QCOMPARE(scanTimeBits(actual.scanTimeFromFrameIndex(frame)),
                     scanTimeBits(expected.scanTimeFromFrameIndex(frame)));
        }
        for (auto it = inputs.times.begin(); it != inputs.times.end(); ++it) {
            for (std::int64_t number : {std::int64_t(it.key()) - 1, std::int64_t(it.key()),
                                       std::int64_t(it.key()) + 1}) {
                if (number < std::numeric_limits<int>::min() || number > std::numeric_limits<int>::max())
                    continue;
                QCOMPARE(scanTimeBits(actual.scanTimeFromScanNumber(int(number))),
                         scanTimeBits(expected.scanTimeFromScanNumber(int(number))));
            }
            for (float query : {it.value(), std::nextafter(it.value(), -std::numeric_limits<float>::infinity()),
                                std::nextafter(it.value(), std::numeric_limits<float>::infinity()),
                                it.value() - .5f, it.value() + .5f}) {
                FrameIndex frame = -999;
                QCOMPARE(actual.frameIndexFromScanTime(query, &frame), eNoError);
                QCOMPARE(frame, expected.frameIndexFromScanTime(query));
            }
        }
        for (float query : {-1e30f, 1e30f, .5f / 60.f, -0.f, 0.f,
                            std::numeric_limits<float>::infinity(),
                            -std::numeric_limits<float>::infinity(),
                            std::numeric_limits<float>::quiet_NaN()}) {
            FrameIndex frame = -999;
            QCOMPARE(actual.frameIndexFromScanTime(query, &frame), eNoError);
            QCOMPARE(frame, expected.frameIndexFromScanTime(query));
        }
    }
}

QTEST_MAIN(MsFrameTests)

#include "MsFrameTests.moc"
