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

    void preprocessMsFrame();


};


void MsFrameTests::preprocessMsFrame() {
    QSKIP("TODO: enable with internal test data");

    ERR_INIT

    QSKIP("activate when proper pathing is used");

    const QString prqfilepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

    MsReaderParquet msReader;
    e = msReader.openFile(prqfilepath);
    QCOMPARE(e, Error::eNoError);

    QMap<ScanNumber, ScanPoints> ms1Scans;

    const int msLevel = 1;
    e = msReader.getScanPoints(
            msLevel,
            &ms1Scans
    );
    QCOMPARE(e, eNoError);

//    MsFrame msFrame;
//    e = msFrame.init(
//            ms1Scans,
//            msReader.getScanNumberVsScanTime()
//            );
//    QCOMPARE(e, eNoError);

//    const bool denoise = true;
//    const bool deisotope = true;
//    const bool smooth = false; //TODO turn on when implemented
//    e = msFrame.preprocessMsFrame(
//            denoise,
//            deisotope,
//            smooth
//            );
//    QCOMPARE(e, eNoError);
//
//    const QMap<ScanNumber, ScanPoints> &ms1ScansProcessed = msFrame.m_frame;
//
//    const QVector<int> expectedScanSizes = {
//        22, 135, 128, 122, 134, 111, 106, 116, 110, 118, 104, 111, 106, 106, 114, 117, 113, 122, 118, 119, 118, 121, 129, 119, 115,
//        115, 121, 118, 118, 123, 111, 122, 104, 108, 127, 110, 113, 110, 123, 108, 100, 108, 112, 98, 103, 103, 110, 99, 113, 109,
//        110, 106, 106, 115, 107, 96, 113, 102, 107, 107, 97, 104, 121, 106, 118, 107, 117, 113, 109, 111, 111, 115, 112, 118, 120,
//        122, 148, 145, 134, 143, 147, 149, 151, 137, 129, 120, 104, 124, 170, 176, 103, 118, 85, 53, 55, 82, 52, 74, 115, 112, 110,
//        112, 129, 96, 103, 107, 140, 161, 196, 145, 156, 114, 142, 126, 121, 116, 122, 129, 167, 180, 158, 164, 125, 131, 143, 169,
//        143, 205, 205, 201, 170, 166, 146, 131, 135, 130, 111, 131, 138, 131, 142, 151, 167, 193, 161, 196, 181, 174, 148, 115, 139,
//        135, 132, 160, 142, 148, 67, 92, 152, 159, 175, 182, 171, 187, 200, 160, 154, 156, 134, 179, 197, 222, 234, 233, 240, 246,
//        219, 177, 149, 115, 144, 188, 184, 204, 217, 230, 222, 210, 190, 194, 196, 221, 216, 215, 253, 218, 209, 192, 225, 247, 228,
//        196, 184, 182, 224, 221, 223, 182, 220, 222, 181, 194, 187, 194, 233, 283, 279, 311, 229, 211, 154, 166, 170, 169, 220, 232,
//        309, 286, 317, 342, 327, 289, 173, 172, 157, 146, 146, 147, 167, 202, 258, 283, 279, 257, 233, 261, 254, 252, 141, 161, 163,
//        197, 260, 292, 329, 301, 272, 297, 276, 254, 269, 213, 249, 259, 274, 285, 342, 314, 242, 253, 213, 277, 223, 225, 230, 230,
//        128, 146, 199, 229, 286, 348, 341, 306, 341, 341, 351, 339, 340, 233, 259, 236, 280, 319, 343, 364, 239, 243, 299, 315, 307,
//        225, 206, 220, 195, 173, 206, 190, 189, 207, 223, 292, 343, 303, 320, 337, 271, 299, 270, 214, 131, 142, 176, 193, 285, 341,
//        359, 347, 308, 189, 227, 271, 303, 262, 281, 303, 244, 236, 76, 79, 127, 204, 278, 279, 308, 307, 245, 244, 199, 108, 166, 249,
//        131, 206, 335, 331, 389, 397, 371, 192, 253, 340, 221, 201, 210, 278, 305, 357, 330, 337, 296, 313, 210, 241, 340, 352, 362, 359,
//        443, 453, 472, 490, 586, 592, 432, 567, 560, 533, 610, 586, 647, 681, 683, 743, 607, 595, 640, 626, 587, 385, 203, 336, 246, 226,
//        224, 212, 140, 248, 596, 648, 684, 630, 493
//    };
//
//    int counter = 0;
//    for (const ScanNumber sn : ms1ScansProcessed.keys()) {
//        QCOMPARE(ms1ScansProcessed.value(sn).size(), expectedScanSizes.at(counter++));
//        QVERIFY(ms1Scans.value(sn).size() > ms1ScansProcessed.value(sn).size());
//    }

}


QTEST_MAIN(MsFrameTests)

#include "MsFrameTests.moc"
