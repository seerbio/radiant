#include <QtTest>
#include <QCoreApplication>

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "EigenUtils.h"
#include "MS2Ion.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"

#include <QVector>

class CandidateScorertronTests : public QObject
{
    Q_OBJECT

public:
    CandidateScorertronTests() = default;
    ~CandidateScorertronTests() override = default;

private slots:
    static void initTest();
    static void integrationDev();
    static void calculateScoresAndOtherStuffTooTest();

private:

    static Err buildInputData(
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanPointsMS1,
            QMap<ScanNumber, ScanPoints*> *ms1FramePtrs,
            MsFrame *msFrameMS1,
            TurboXIC *turboXICMS1
    ) {

        ERR_INIT

        const QString prqFFFilePath
                = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

        e = msReaderPointerAcc->openFile(prqFFFilePath); ree;

        msReaderPointerAcc->ptr->collateMS2MzTargetFrames(diaTargetFrames);

        const int msLevel = 1;
        e = msReaderPointerAcc->ptr->getScanPoints(msLevel, scanNumberVsScanPointsMS1); ree;

        for (auto it = scanNumberVsScanPointsMS1->begin(); it != scanNumberVsScanPointsMS1->end(); it++) {
            ms1FramePtrs->insert(it.key(), &it.value());
        }

        e = msFrameMS1->init(*ms1FramePtrs, msReaderPointerAcc->ptr->getScanNumberVsScanTime()); ree;
        e = turboXICMS1->init(msFrameMS1->frameIndexVsScanPoints()); ree;

        ERR_RETURN
    }

    static TargetDecoyCandidatePair buildTargetDecoyCandidatePair() {

        MS2Ion i1;
        i1.mz = 506.335;
        i1.intensity = 1.0;
        i1.charge = 1;
        i1.ionLabel = "y4";

        MS2Ion i2;
        i2.mz = 407.266;
        i2.intensity = 0.488;
        i2.charge = 1;
        i2.ionLabel = "y3";

        MS2Ion i3;
        i3.mz = 272.125;
        i3.intensity = 0.30;
        i3.charge = 1;
        i3.ionLabel = "b4";

        QVector<MS2Ion> targetIons;

        TargetDecoyCandidatePair targetDecoyCandidatePair(
                PeptideStringWithMods("DGVVLFK"),
                {i1, i2, i3},
                {i1, i2, i3},
                2,
                776.443,
                69.0608,
                3
                );

        return targetDecoyCandidatePair;
    }

    static Err buildFrame(QMap<ScanNumber, ScanPoints> *frame) {

        ERR_INIT

        const QString framePath
                = QDir(qApp->applicationDirPath()).filePath("framePoints.frame");

        QVector<MsFrameScanPointRows> rows;
        e = ParquetReader::read(framePath, &rows); ree;

        for (const MsFrameScanPointRows &r : rows) {

            ScanPoints scanPoints;
            for (int i = 0; i < r.mzVals.size(); i++) {
                ScanPoint scanPoint(r.mzVals.at(i), r.intensityVals.at(i));
                scanPoints.push_back(scanPoint);
            }

            frame->insert(r.frameIndex, scanPoints);
        }

        ERR_RETURN
    }

};


namespace {

    QVector<QPair<float, float>> prodVec = {
        QPair( 0.0155471 , 0 ), QPair( 0.0408375 , 0 ), QPair( 0.0660899 , 0 ), QPair( 0.0918517 , 0 ),
        QPair( 0.117471 , 0 ), QPair( 0.142704 , 0 ),  QPair( 0.168549 , 0 ), QPair( 0.193802 , 0 ), QPair( 0.21955 , 0 ),
        QPair( 0.24488 , 0 ), QPair( 0.270689 , 0 ),    QPair( 0.295972 , 0 ),    QPair( 0.32162 , 0 ),    QPair( 0.346769 , 0 ),    QPair( 0.372569 , 0 ),    QPair( 0.39791 , 0 ),    QPair( 0.423669 , 0 ),    QPair( 0.448999 , 0 ),    QPair( 0.474803 , 0 ),    QPair( 0.500229 , 0 ),    QPair( 0.526018 , 0 ),    QPair( 0.551387 , 0 ),    QPair( 0.577085 , 0 ),    QPair( 0.602357 , 0 ),    QPair( 0.628134 , 0 ),    QPair( 0.653508 , 0 ),    QPair( 0.67935 , 0 ),    QPair( 0.704672 , 0 ),    QPair( 0.730548 , 0 ),    QPair( 0.755991 , 0 ),    QPair( 0.78189 , 0 ),    QPair( 0.80738 , 0 ),    QPair( 0.833103 , 0 ),    QPair( 0.858458 , 0 ),    QPair( 0.884358 , 0 ),    QPair( 0.90964 , 0 ),    QPair( 0.935523 , 0 ),    QPair( 0.960932 , 0 ),    QPair( 0.986673 , 0 ),    QPair( 1.01207 , 0 ),    QPair( 1.03779 , 0 ),    QPair( 1.06311 , 0 ),    QPair( 1.08885 , 0 ),    QPair( 1.11418 , 0 ),    QPair( 1.13999 , 0 ),    QPair( 1.16515 , 0 ),    QPair( 1.19103 , 0 ),    QPair( 1.2162 , 0 ),    QPair( 1.24197 , 0 ),    QPair( 1.26726 , 0 ),    QPair( 1.29308 , 0 ),    QPair( 1.31845 , 0 ),    QPair( 1.34425 , 0 ),    QPair( 1.37006 , 0 ),    QPair( 1.39537 , 0 ),    QPair( 1.42122 , 0 ),    QPair( 1.44648 , 0 ),    QPair( 1.47228 , 0 ),    QPair( 1.49756 , 0 ),    QPair( 1.52335 , 0 ),    QPair( 1.54869 , 0 ),    QPair( 1.57431 , 0 ),    QPair( 1.59962 , 0 ),    QPair( 1.62531 , 0 ),    QPair( 1.65057 , 0 ),    QPair( 1.6763 , 0 ),    QPair( 1.70159 , 0 ),    QPair( 1.72742 , 0 ),    QPair( 1.75279 , 0 ),    QPair( 1.77861 , 0 ),    QPair( 1.80399 , 0 ),    QPair( 1.82991 , 0 ),    QPair( 1.85529 , 0 ),    QPair( 1.88111 , 0 ),    QPair( 1.90641 , 0 ),    QPair( 1.93215 , 0 ),    QPair( 1.95755 , 0 ),    QPair( 1.98339 , 0 ),    QPair( 2.00868 , 0 ),    QPair( 2.03446 , 0 ),    QPair( 2.05972 , 0 ),    QPair( 2.08556 , 0 ),    QPair( 2.11085 , 0 ),    QPair( 2.13661 , 0 ),    QPair( 2.16187 , 0 ),    QPair( 2.18748 , 0 ),    QPair( 2.21273 , 0 ),    QPair( 2.23849 , 0 ),    QPair( 2.26375 , 0 ),    QPair( 2.28948 , 0 ),    QPair( 2.31494 , 0 ),    QPair( 2.34062 , 0 ),    QPair( 2.36581 , 0 ),    QPair( 2.39158 , 0 ),    QPair( 2.41693 , 0 ),    QPair( 2.44273 , 0 ),    QPair( 2.46799 , 0 ),    QPair( 2.49383 , 0 ),    QPair( 2.5191 , 0 ),    QPair( 2.54493 , 0 ),    QPair( 2.57031 , 0 ),    QPair( 2.59609 , 0 ),    QPair( 2.62119 , 0 ),    QPair( 2.64685 , 0 ),    QPair( 2.67218 , 0 ),    QPair( 2.69806 , 0 ),    QPair( 2.72332 , 0 ),    QPair( 2.74903 , 0 ),    QPair( 2.77426 , 0 ),    QPair( 2.8 , 0 ),    QPair( 2.82523 , 0 ),    QPair( 2.85097 , 0 ),    QPair( 2.87634 , 0 ),    QPair( 2.90211 , 0 ),    QPair( 2.92745 , 0 ),    QPair( 2.95323 , 0 ),    QPair( 2.97909 , 0 ),    QPair( 3.00438 , 0 ),    QPair( 3.03021 , 0 ),    QPair( 3.05537 , 0 ),    QPair( 3.08123 , 0 ),    QPair( 3.10645 , 0 ),    QPair( 3.13208 , 0 ),    QPair( 3.15742 , 0 ),    QPair( 3.18332 , 0 ),    QPair( 3.20863 , 0 ),    QPair( 3.23437 , 0 ),    QPair( 3.26175 , 0 ),    QPair( 3.28753 , 0 ),    QPair( 3.31287 , 0 ),    QPair( 3.33843 , 0 ),    QPair( 3.36378 , 0 ),    QPair( 3.38947 , 0 ),    QPair( 3.4154 , 0 ),    QPair( 3.44081 , 0 ),    QPair( 3.46671 , 0 ),    QPair( 3.49193 , 0 ),    QPair( 3.51778 , 0 ),    QPair( 3.54333 , 0 ),    QPair( 3.56935 , 0 ),    QPair( 3.59485 , 0 ),    QPair( 3.62075 , 0 ),    QPair( 3.64626 , 0 ),    QPair( 3.67218 , 0 ),    QPair( 3.69771 , 0 ),    QPair( 3.72349 , 0 ),    QPair( 3.74883 , 0 ),    QPair( 3.7747 , 0 ),    QPair( 3.8 , 0 ),    QPair( 3.82574 , 0 ),    QPair( 3.85078 , 0 ),    QPair( 3.87642 , 0 ),    QPair( 3.90162 , 0 ),    QPair( 3.92711 , 0 ),    QPair( 3.95209 , 0 ),    QPair( 3.97748 , 0 ),    QPair( 4.00253 , 0 ),    QPair( 4.0281 , 0 ),    QPair( 4.05324 , 0 ),    QPair( 4.07878 , 0 ),    QPair( 4.10395 , 0 ),    QPair( 4.12969 , 0 ),    QPair( 4.1548 , 0 ),    QPair( 4.18041 , 0 ),    QPair( 4.20551 , 0 ),    QPair( 4.23116 , 0 ),    QPair( 4.25618 , 0 ),    QPair( 4.28175 , 0 ),    QPair( 4.30699 , 0 ),    QPair( 4.33236 , 0 ),    QPair( 4.3575 , 0 ),    QPair( 4.38295 , 0 ),    QPair( 4.40806 , 0 ),    QPair( 4.43336 , 0 ),    QPair( 4.45848 , 0 ),    QPair( 4.484 , 0 ),    QPair( 4.50927 , 0 ),    QPair( 4.5346 , 0 ),    QPair( 4.55944 , 0 ),    QPair( 4.58486 , 0 ),    QPair( 4.60963 , 0 ),    QPair( 4.63435 , 0 ),    QPair( 4.65937 , 0 ),    QPair( 4.68386 , 0 ),    QPair( 4.7086 , 0 ),    QPair( 4.73303 , 0 ),    QPair( 4.75787 , 0 ),    QPair( 4.78231 , 0 ),    QPair( 4.80672 , 0 ),    QPair( 4.83163 , 0 ),    QPair( 4.85594 , 0 ),    QPair( 4.88078 , 0 ),    QPair( 4.90522 , 0 ),    QPair( 4.9299 , 0 ),    QPair( 4.95507 , 0 ),    QPair( 4.97955 , 0 ),    QPair( 5.00452 , 0 ),    QPair( 5.02919 , 0 ),    QPair( 5.05429 , 0 ),    QPair( 5.07903 , 0 ),    QPair( 5.10406 , 0 ),    QPair( 5.12879 , 0 ),    QPair( 5.15402 , 0 ),    QPair( 5.17854 , 0 ),    QPair( 5.20291 , 0 ),    QPair( 5.2275 , 0 ),    QPair( 5.25193 , 0 ),    QPair( 5.27682 , 0 ),    QPair( 5.3013 , 0 ),    QPair( 5.32585 , 0 ),    QPair( 5.35086 , 0 ),    QPair( 5.37514 , 0 ),    QPair( 5.39993 , 0 ),    QPair( 5.4246 , 0 ),    QPair( 5.44961 , 0 ),    QPair( 5.47398 , 0 ),    QPair( 5.49837 , 0 ),    QPair( 5.52303 , 0 ),    QPair( 5.5474 , 0 ),    QPair( 5.57211 , 0 ),    QPair( 5.59671 , 0 ),    QPair( 5.62115 , 0 ),    QPair( 5.64634 , 0 ),    QPair( 5.67083 , 0 ),    QPair( 5.69596 , 0 ),    QPair( 5.72057 , 0 ),    QPair( 5.74567 , 0 ),    QPair( 5.77027 , 0 ),    QPair( 5.79529 , 0 ),    QPair( 5.81995 , 0 ),    QPair( 5.84463 , 0 ),    QPair( 5.86983 , 0 ),    QPair( 5.89429 , 0 ),    QPair( 5.91952 , 0 ),    QPair( 5.94431 , 0 ),    QPair( 5.96954 , 0 ),    QPair( 5.99439 , 0 ),    QPair( 6.01955 , 0 ),    QPair( 6.04431 , 0 ),    QPair( 6.06918 , 0 ),    QPair( 6.09444 , 0 ),    QPair( 6.11934 , 0 ),    QPair( 6.14457 , 0 ),    QPair( 6.16958 , 0 ),    QPair( 6.19488 , 0 ),    QPair( 6.2196 , 0 ),    QPair( 6.24478 , 0 ),    QPair( 6.26955 , 0 ),    QPair( 6.29478 , 0 ),    QPair( 6.31946 , 0 ),    QPair( 6.34451 , 0 ),    QPair( 6.36929 , 0 ),    QPair( 6.39384 , 0 ),    QPair( 6.41852 , 0 ),    QPair( 6.44289 , 0 ),    QPair( 6.46775 , 0 ),    QPair( 6.49203 , 0 ),    QPair( 6.51653 , 0 ),    QPair( 6.54139 , 0 ),    QPair( 6.56633 , 0 ),    QPair( 6.59147 , 0 ),    QPair( 6.61621 , 0 ),    QPair( 6.64142 , 0 ),    QPair( 6.66603 , 0 ),    QPair( 6.69115 , 0 ),    QPair( 6.71579 , 0 ),    QPair( 6.74048 , 0 ),    QPair( 6.76544 , 0 ),    QPair( 6.79009 , 0 ),    QPair( 6.81515 , 0 ),    QPair( 6.83993 , 0 ),    QPair( 6.86506 , 0 ),    QPair( 6.88984 , 0 ),    QPair( 6.91467 , 0 ),    QPair( 6.93979 , 0 ),    QPair( 6.96467 , 0 ),    QPair( 6.98988 , 0 ),    QPair( 7.01441 , 0 ),    QPair( 7.03952 , 0 ),    QPair( 7.06439 , 0 ),    QPair( 7.08942 , 0 ),    QPair( 7.11412 , 0 ),    QPair( 7.13881 , 0 ),    QPair( 7.16389 , 0 ),    QPair( 7.18853 , 0 ),    QPair( 7.21375 , 0 ),    QPair( 7.23848 , 0 ),    QPair( 7.2637 , 0 ),    QPair( 7.28834 , 0 ),    QPair( 7.31279 , 0 ),    QPair( 7.33767 , 0 ),    QPair( 7.36223 , 0 ),    QPair( 7.38729 , 0 ),    QPair( 7.41202 , 0 ),    QPair( 7.43705 , 0 ),    QPair( 7.46162 , 0 ),    QPair( 7.48669 , 0 ),    QPair( 7.51136 , 0 ),    QPair( 7.53614 , 0 ),    QPair( 7.56113 , 0 ),    QPair( 7.58563 , 0 ),    QPair( 7.61061 , 0 ),    QPair( 7.63531 , 0 ),    QPair( 7.66045 , 0 ),    QPair( 7.68524 , 0 ),    QPair( 7.7103 , 0 ),    QPair( 7.73481 , 0 ),    QPair( 7.75935 , 0 ),    QPair( 7.78433 , 0 ),    QPair( 7.80897 , 0 ),    QPair( 7.83399 , 0 ),    QPair( 7.85851 , 0 ),    QPair( 7.88351 , 0 ),    QPair( 7.90785 , 0 ),    QPair( 7.93226 , 0 ),    QPair( 7.95714 , 0 ),    QPair( 7.98157 , 0 ),    QPair( 8.00633 , 0 ),    QPair( 8.03081 , 0 ),    QPair( 8.05586 , 0 ),    QPair( 8.08061 , 0 ),    QPair( 8.10586 , 0 ),    QPair( 8.13069 , 0 ),    QPair( 8.15555 , 0 ),    QPair( 8.18068 , 0 ),    QPair( 8.20493 , 0 ),    QPair( 8.23186 , 0 ),    QPair( 8.25587 , 0 ),    QPair( 8.28037 , 0 ),    QPair( 8.30486 , 0 ),    QPair( 8.3294 , 0 ),    QPair( 8.35415 , 0 ),    QPair( 8.37866 , 0 ),    QPair( 8.40342 , 0 ),    QPair( 8.42789 , 0 ),    QPair( 8.45277 , 0 ),    QPair( 8.47714 , 0 ),    QPair( 8.50175 , 0 ),    QPair( 8.52682 , 0 ),    QPair( 8.5515 , 0 ),    QPair( 8.57656 , 0 ),    QPair( 8.60139 , 0 ),    QPair( 8.62638 , 0 ),    QPair( 8.65104 , 0 ),    QPair( 8.67606 , 0 ),    QPair( 8.70061 , 0 ),    QPair( 8.72518 , 0 ),    QPair( 8.75047 , 0 ),    QPair( 8.775 , 0 ),    QPair( 8.80025 , 0 ),    QPair( 8.82503 , 0 ),    QPair( 8.85012 , 0 ),    QPair( 8.87473 , 0 ),    QPair( 8.89931 , 0 ),    QPair( 8.92433 , 0 ),    QPair( 8.94898 , 0 ),    QPair( 8.9739 , 0 ),    QPair( 8.99836 , 0 ),    QPair( 9.02322 , 0 ),    QPair( 9.04755 , 0 ),    QPair( 9.07188 , 0 ),    QPair( 9.09684 , 0 ),    QPair( 9.12145 , 0 ),    QPair( 9.14648 , 0 ),    QPair( 9.1712 , 0 ),    QPair( 9.19623 , 0 ),    QPair( 9.221 , 0 ),    QPair( 9.24569 , 0 ),    QPair( 9.27077 , 0 ),    QPair( 9.29543 , 0 ),    QPair( 9.32037 , 0 ),    QPair( 9.34475 , 0 ),    QPair( 9.36964 , 0 ),    QPair( 9.39429 , 0 ),    QPair( 9.41924 , 0 ),    QPair( 9.44378 , 0 ),    QPair( 9.46832 , 0 ),    QPair( 9.49326 , 0 ),    QPair( 9.51796 , 0 ),    QPair( 9.5429 , 0 ),    QPair( 9.5674 , 0 ),    QPair( 9.59206 , 0 ),    QPair( 9.61646 , 0 ),    QPair( 9.64082 , 0 ),    QPair( 9.66556 , 0 ),    QPair( 9.69003 , 0 ),    QPair( 9.71476 , 0 ),    QPair( 9.73917 , 0 ),    QPair( 9.76371 , 0 ),    QPair( 9.78862 , 0 ),    QPair( 9.81321 , 0 ),    QPair( 9.83826 , 0 ),    QPair( 9.86296 , 0 ),    QPair( 9.88807 , 0 ),    QPair( 9.91258 , 0 ),    QPair( 9.93748 , 0 ),    QPair( 9.96198 , 0 ),    QPair( 9.98685 , 0 ),    QPair( 10.012 , 0 ),    QPair( 10.0368 , 0 ),    QPair( 10.0619 , 0 ),    QPair( 10.0867 , 0 ),    QPair( 10.1118 , 0 ),    QPair( 10.1366 , 0 ),    QPair( 10.1618 , 0 ),    QPair( 10.1864 , 0 ),    QPair( 10.2111 , 0 ),    QPair( 10.2362 , 0 ),    QPair( 10.2605 , 0 ),    QPair( 10.2856 , 0 ),    QPair( 10.3104 , 0 ),    QPair( 10.3356 , 0 ),    QPair( 10.3603 , 0 ),    QPair( 10.3854 , 0 ),    QPair( 10.4103 , 0 ),    QPair( 10.435 , 0 ),    QPair( 10.46 , 0 ),    QPair( 10.4847 , 0 ),    QPair( 10.5098 , 0 ),    QPair( 10.5344 , 0 ),    QPair( 10.5595 , 0 ),    QPair( 10.5838 , 0 ),    QPair( 10.6079 , 0 ),    QPair( 10.6327 , 0 ),    QPair( 10.657 , 0 ),    QPair( 10.6816 , 0 ),    QPair( 10.7061 , 0 ),    QPair( 10.7306 , 0 ),    QPair( 10.7553 , 0 ),    QPair( 10.7798 , 0 ),    QPair( 10.8049 , 0 ),    QPair( 10.8296 , 0 ),    QPair( 10.8547 , 0 ),    QPair( 10.8796 , 0 ),    QPair( 10.9044 , 0 ),    QPair( 10.9299 , 0 ),    QPair( 10.9549 , 0 ),    QPair( 10.9802 , 0 ),    QPair( 11.0049 , 0 ),    QPair( 11.03 , 0 ),    QPair( 11.0546 , 0 ),    QPair( 11.0797 , 0 ),    QPair( 11.1044 , 0 ),    QPair( 11.1298 , 0 ),    QPair( 11.1546 , 0 ),    QPair( 11.1796 , 0 ),    QPair( 11.2042 , 0 ),    QPair( 11.2289 , 0 ),    QPair( 11.2538 , 0 ),    QPair( 11.2781 , 0 ),    QPair( 11.303 , 0 ),    QPair( 11.3277 , 0 ),    QPair( 11.3528 , 0 ),    QPair( 11.3772 , 0 ),    QPair( 11.4015 , 0 ),    QPair( 11.4263 , 0 ),    QPair( 11.4507 , 0 ),    QPair( 11.4756 , 0 ),    QPair( 11.5001 , 0 ),    QPair( 11.5246 , 0 ),    QPair( 11.5495 , 0 ),    QPair( 11.5741 , 0 ),    QPair( 11.5992 , 0 ),    QPair( 11.6239 , 0 ),    QPair( 11.6491 , 0 ),    QPair( 11.6738 , 0 ),    QPair( 11.6991 , 0 ),    QPair( 11.724 , 0 ),    QPair( 11.7491 , 0 ),    QPair( 11.7737 , 0 ),    QPair( 11.7982 , 0 ),    QPair( 11.823 , 0 ),    QPair( 11.8476 , 0 ),    QPair( 11.8725 , 0 ),    QPair( 11.897 , 0 ),    QPair( 11.9214 , 0 ),    QPair( 11.946 , 0 ),    QPair( 11.9705 , 0 ),    QPair( 11.9954 , 0 ),    QPair( 12.0199 , 0 ),    QPair( 12.045 , 0 ),    QPair( 12.0699 , 0 ),    QPair( 12.0946 , 0 ),    QPair( 12.1197 , 0 ),    QPair( 12.1445 , 0 ),    QPair( 12.1698 , 0 ),    QPair( 12.1947 , 0 ),    QPair( 12.2199 , 0 ),    QPair( 12.2448 , 0 ),    QPair( 12.2699 , 0 ),    QPair( 12.2947 , 0 ),    QPair( 12.3198 , 0 ),    QPair( 12.3447 , 0 ),    QPair( 12.3697 , 0 ),    QPair( 12.3947 , 0 ),    QPair( 12.4195 , 0 ),    QPair( 12.4447 , 0 ),    QPair( 12.4694 , 0 ),    QPair( 12.4945 , 0 ),    QPair( 12.5192 , 0 ),    QPair( 12.5443 , 0 ),    QPair( 12.5689 , 0 ),    QPair( 12.5938 , 0 ),    QPair( 12.6183 , 0 ),    QPair( 12.6427 , 0 ),    QPair( 12.6677 , 0 ),    QPair( 12.6923 , 0 ),    QPair( 12.7174 , 0 ),    QPair( 12.7421 , 0 ),    QPair( 12.7673 , 0 ),    QPair( 12.7921 , 0 ),    QPair( 12.8172 , 0 ),    QPair( 12.8419 , 0 ),    QPair( 12.8669 , 0 ),    QPair( 12.8917 , 0 ),    QPair( 12.9165 , 0 ),    QPair( 12.9416 , 0 ),    QPair( 12.9663 , 0 ),    QPair( 12.9915 , 0 ),    QPair( 13.0164 , 0 ),    QPair( 13.0416 , 0 ),    QPair( 13.0664 , 0 ),    QPair( 13.0915 , 0 ),    QPair( 13.1161 , 0 ),    QPair( 13.1413 , 0 ),    QPair( 13.1661 , 0 ),    QPair( 13.1933 , 0 ),    QPair( 13.2182 , 0 ),    QPair( 13.243 , 0 ),    QPair( 13.2682 , 0 ),    QPair( 13.2929 , 0 ),    QPair( 13.3183 , 0 ),    QPair( 13.3428 , 0 ),    QPair( 13.3679 , 0 ),    QPair( 13.3924 , 0 ),    QPair( 13.4174 , 0 ),    QPair( 13.442 , 0 ),    QPair( 13.467 , 0 ),    QPair( 13.4923 , 0 ),    QPair( 13.5169 , 0 ),    QPair( 13.542 , 0 ),    QPair( 13.5665 , 0 ),    QPair( 13.5917 , 0 ),    QPair( 13.6163 , 0 ),    QPair( 13.641 , 0 ),    QPair( 13.666 , 0 ),    QPair( 13.6908 , 0 ),    QPair( 13.716 , 0 ),    QPair( 13.7408 , 0 ),    QPair( 13.766 , 0 ),    QPair( 13.7908 , 0 ),    QPair( 13.816 , 0 ),    QPair( 13.8408 , 0 ),    QPair( 13.8661 , 0 ),    QPair( 13.8911 , 0 ),    QPair( 13.9161 , 0 ),    QPair( 13.9408 , 0 ),    QPair( 13.9656 , 0 ),    QPair( 13.9909 , 0 ),    QPair( 14.0156 , 0 ),    QPair( 14.0405 , 0 ),    QPair( 14.0652 , 0 ),    QPair( 14.0904 , 0 ),    QPair( 14.1152 , 0 ),    QPair( 14.1402 , 0 ),    QPair( 14.165 , 0 ),    QPair( 14.1905 , 0 ),    QPair( 14.2154 , 0 ),    QPair( 14.2404 , 0 ),    QPair( 14.2658 , 0 ),    QPair( 14.2907 , 0 ),    QPair( 14.3162 , 0.437702 ),    QPair( 14.3413 , 1.17127 ),    QPair( 14.3666 , 2.19535 ),    QPair( 14.3916 , 3.1629 ),    QPair( 14.4169 , 3.8947 ),    QPair( 14.4417 , 3.14549 ),    QPair( 14.467 , 2.4013 ),    QPair( 14.4918 , 1.66397 ),    QPair( 14.5165 , 0.7364 ),    QPair( 14.5414 , 0 ),    QPair( 14.5659 , 0 ),    QPair( 14.591 , 0 ),    QPair( 14.6155 , 0 ),    QPair( 14.6405 , 0 ),    QPair( 14.6649 , 0 ),    QPair( 14.6898 , 0 ),    QPair( 14.7145 , 0.741716 ),    QPair( 14.7392 , 1.48959 ),    QPair( 14.7643 , 2.23438 ),    QPair( 14.7892 , 2.23896 ),    QPair( 14.8144 , 2.25047 ),    QPair( 14.8391 , 1.49399 ),    QPair( 14.8641 , 0.745565 ),    QPair( 14.8887 , 0 ),    QPair( 14.9133 , 0.743481 ),    QPair( 14.9384 , 1.68433 ),    QPair( 14.963 , 2.63828 ),    QPair( 14.9883 , 3.82186 ),    QPair( 15.0132 , 5.02671 ),    QPair( 15.0383 , 5.4313 ),    QPair( 15.063 , 5.6025 ),    QPair( 15.088 , 5.79021 ),    QPair( 15.1127 , 5.74244 ),    QPair( 15.138 , 5.7227 ),    QPair( 15.163 , 5.69258 ),    QPair( 15.1879 , 5.67803 ),    QPair( 15.2133 , 5.4513 ),    QPair( 15.2384 , 5.23616 ),    QPair( 15.2636 , 4.85579 ),    QPair( 15.2886 , 4.89248 ),    QPair( 15.3139 , 4.55568 ),    QPair( 15.3387 , 4.55882 ),    QPair( 15.3639 , 4.39057 ),    QPair( 15.3887 , 3.62074 ),    QPair( 15.4138 , 3.60085 ),    QPair( 15.4386 , 3.74588 ),    QPair( 15.4639 , 3.94905 ),    QPair( 15.4887 , 4.3336 ),    QPair( 15.5141 , 5.50352 ),    QPair( 15.5391 , 5.50556 ),    QPair( 15.564 , 5.7171 ),    QPair( 15.5893 , 5.71515 ),    QPair( 15.6143 , 5.71223 ),    QPair( 15.6396 , 5.67928 ),    QPair( 15.6644 , 5.66766 ),    QPair( 15.6896 , 5.66564 ),    QPair( 15.7145 , 5.66438 ),    QPair( 15.7398 , 5.65671 ),    QPair( 15.7648 , 5.64943 ),    QPair( 15.79 , 5.64268 ),    QPair( 15.8149 , 5.64285 ),    QPair( 15.8402 , 5.64056 ),    QPair( 15.8647 , 5.64366 ),    QPair( 15.8891 , 5.64458 ),    QPair( 15.914 , 5.64774 ),    QPair( 15.9387 , 5.6461 ),    QPair( 15.9638 , 5.64483 ),    QPair( 15.988 , 5.64188 ),    QPair( 16.0126 , 5.64106 ),    QPair( 16.037 , 5.64514 ),    QPair( 16.0614 , 5.63694 ),    QPair( 16.0863 , 5.62664 ),    QPair( 16.1107 , 5.61279 ),    QPair( 16.1356 , 5.5968 ),    QPair( 16.1602 , 5.57654 ),    QPair( 16.1847 , 5.57558 ),    QPair( 16.2096 , 5.57072 ),    QPair( 16.234 , 5.5653 ),    QPair( 16.259 , 5.56385 ),    QPair( 16.2836 , 5.57264 ),    QPair( 16.3088 , 5.58234 ),    QPair( 16.3333 , 5.59628 ),    QPair( 16.3579 , 5.61677 ),    QPair( 16.3828 , 5.63912 ),    QPair( 16.4074 , 5.66625 ),    QPair( 16.4323 , 5.69825 ),    QPair( 16.4568 , 5.73271 ),    QPair( 16.4814 , 5.76173 ),    QPair( 16.5063 , 5.78425 ),    QPair( 16.531 , 5.79006 ),    QPair( 16.5559 , 5.78491 ),    QPair( 16.5805 , 5.77646 ),    QPair( 16.6055 , 5.76947 ),    QPair( 16.6301 , 5.75703 ),    QPair( 16.6551 , 5.73783 ),    QPair( 16.6799 , 5.71869 ),    QPair( 16.7048 , 5.7136 ),    QPair( 16.7302 , 5.71108 ),    QPair( 16.7551 , 5.7074 ),    QPair( 16.7805 , 5.69314 ),    QPair( 16.8056 , 5.66192 ),    QPair( 16.831 , 5.61362 ),    QPair( 16.8559 , 5.53224 ),    QPair( 16.8811 , 5.47738 ),    QPair( 16.9058 , 5.40783 ),    QPair( 16.931 , 5.38065 ),    QPair( 16.9557 , 5.3719 ),    QPair( 16.9804 , 5.4035 ),    QPair( 17.0054 , 5.43646 ),    QPair( 17.0303 , 5.51329 ),    QPair( 17.0554 , 5.40068 ),    QPair( 17.0802 , 5.45168 ),    QPair( 17.1056 , 5.47503 ),    QPair( 17.1306 , 5.27545 ),    QPair( 17.156 , 5.23898 ),    QPair( 17.1809 , 5.38046 ),    QPair( 17.2061 , 5.33933 ),    QPair( 17.231 , 5.13712 ),    QPair( 17.2557 , 5.13001 ),    QPair( 17.2808 , 4.96005 ),    QPair( 17.3056 , 4.78168 ),    QPair( 17.3309 , 4.60448 ),    QPair( 17.3557 , 4.78953 ),    QPair( 17.3808 , 4.98473 ),    QPair( 17.4056 , 5.18933 ),    QPair( 17.4308 , 5.18203 ),    QPair( 17.4558 , 5.17313 ),    QPair( 17.4811 , 5.1873 ),    QPair( 17.5059 , 5.21468 ),    QPair( 17.5304 , 5.22614 ),    QPair( 17.5556 , 5.4138 ),    QPair( 17.5802 , 5.60532 ),    QPair( 17.6053 , 5.62426 ),    QPair( 17.6299 , 5.64519 ),    QPair( 17.655 , 5.67883 ),    QPair( 17.6798 , 5.74927 ),    QPair( 17.705 , 5.61675 ),    QPair( 17.7299 , 5.43544 ),    QPair( 17.7548 , 5.2325 ),    QPair( 17.7801 , 4.96672 ),    QPair( 17.8048 , 4.69887 ),    QPair( 17.83 , 4.59969 ),    QPair( 17.8548 , 4.54695 ),    QPair( 17.88 , 4.48032 ),    QPair( 17.9048 , 4.46573)
        };

    Err simpleIntegrator(
        const Eigen::VectorX<float> &vec,
        float stopThresholdFraction,
        QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(vec.size() > 0); ree;

        Eigen::VectorX<float> eVec = vec;
        EigenUtils::thresholdVector(static_cast<float>(1.01), &eVec);

        const QMap<int, float> vecApexs = EigenUtils::apexes(eVec);
        if (vecApexs.isEmpty()) {
            ERR_RETURN
        }

        Eigen::VectorX<float> apexes =EigenUtils::convertQMapToEigenVector(vecApexs, vecApexs.lastKey() + 1);
        QVector<QPair<int, float>> apexPairs = EigenUtils::returnTopXIndexAndValues(apexes, vecApexs.size());

        for (const QPair<int, float> &pr : apexPairs) {

            const int apexIndex = pr.first;
            const float apexValue = pr.second;

            if (MathUtils::tZero(apexValue)) {
                continue;
            }

            const float stopThreshold = apexValue * stopThresholdFraction;

            float rightStopVal = apexValue;
            int rightStopIndex = apexIndex;

            int rightCurrentIndex = apexIndex;
            while (rightCurrentIndex < eVec.size()) {

                const float currentValue = eVec(rightCurrentIndex);
                if (currentValue < stopThreshold) {
                    rightStopIndex = rightCurrentIndex;
                    break;
                }

                if (currentValue < rightStopVal || MathUtils::tSame(currentValue, rightStopVal)) {
                    rightStopVal = currentValue;
                    rightStopIndex = rightCurrentIndex;
                    ++rightCurrentIndex;
                    continue;
                }

                break;
            }

            float leftStopVal = apexValue;
            int leftStopIndex = apexIndex;

            int leftCurrentIndex = apexIndex;
            while (leftCurrentIndex < eVec.size()) {

                const float currentValue = eVec(leftCurrentIndex);
                if (currentValue < stopThreshold) {
                    leftStopIndex = leftCurrentIndex;
                    break;
                }

                if (currentValue < leftStopVal || MathUtils::tSame(currentValue, leftStopVal)) {
                    leftStopVal = currentValue;
                    leftStopIndex = leftCurrentIndex;
                    --leftCurrentIndex;
                    continue;
                }

                break;
            }

            peakIntegrationIndexesVsIntensity->push_back({
                 {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size() - 1))},
                 apexValue
            }
            );

            for (int i = leftStopIndex; i <= rightStopIndex; i++) {
                eVec.coeffRef(i) = 0.0;
            }
        }

        std::sort(
            peakIntegrationIndexesVsIntensity->rbegin(),
            peakIntegrationIndexesVsIntensity->rend(),
            [](const QPair<PeakIntegrationIndexes, float> &l, const QPair<PeakIntegrationIndexes, float> &r) {
                return l.second < r.second;
            }
            );

        ERR_RETURN
    }

}//namespace
void CandidateScorertronTests::integrationDev() {

    ERR_INIT

    QVector<float> rt;
    std::transform(
        prodVec.begin(),
        prodVec.end(),
        std::back_inserter(rt),
        [](const QPair<float, float> &pr){return pr.first;}
        );

    QVector<float> prod;
    std::transform(
        prodVec.begin(),
        prodVec.end(),
        std::back_inserter(prod),
        [](const QPair<float, float> &pr){return pr.second;}
        );



    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = simpleIntegrator(
        EigenUtils::convertQVectorToEigenVector(prod),
        0.65f,
        &peakIntegrationsVsIntensities
        );

    for (const auto pr : peakIntegrationsVsIntensities) {
        qDebug() << pr << rt.at(pr.first.first) << rt.at(pr.first.second);
    }

}

void CandidateScorertronTests::initTest() {

    ERR_INIT

    QSKIP("Fix tests for new API");

    // MsReaderPointerAcc msReaderPointerAcc;
    // QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    // QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
    // QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    // MsFrame msFrameMS1;
    // TurboXIC turboXICMS1;
    //
    // e = buildInputData(
    //         &msReaderPointerAcc,
    //         &diaTargetFrames,
    //         &scanNumberVsScanPointsMS1,
    //         &ms1FramePtrs,
    //         &msFrameMS1,
    //         &turboXICMS1
    //         );
    // QCOMPARE(e, eNoError);
    //
    // MsCalibratomatic msCalibratomatic;
    // TurboXIC dummyXIC;
    //
    // CandidateScorertron candidateScorertron;
    // e = candidateScorertron.init(
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         msCalibratomatic,
    //         "100",
    //         8,
    //         3,
    //         &msCalibratomatic,
    //         &dummyXIC
    //         );
    // QCOMPARE(e, eEmptyContainerError);
    //
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameters(),
    //         "",
    //         {},
    //         6,
    //         {},
    //         &msCalibratomatic,
    //         &dummyXIC
    // );
    // QCOMPARE(e, eError);
    //
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         "",
    //         {},
    //         6,
    //         {},
    //         &msCalibratomatic,
    //         &dummyXIC
    // );
    // QCOMPARE(e, eEmptyContainerError);
    //
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         diaTargetFrames.firstKey(),
    //         {},
    //         6,
    //         {},
    //         &msCalibratomatic,
    //         &dummyXIC
    // );
    // QCOMPARE(e, eError);
    //
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         diaTargetFrames.firstKey(),
    //         {10, 10},
    //         -1,
    //         {},
    //         &msCalibratomatic,
    //         &dummyXIC
    // );
    // QCOMPARE(e, eError);
    //
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         diaTargetFrames.firstKey(),
    //         {10, 10},
    //         6,
    //         {},
    //         &msCalibratomatic,
    //         &dummyXIC
    // );
    // QCOMPARE(e, eError);
    //
    // const QHash<MzHashed, int> mzHashedVsCount;
    // e = candidateScorertron.init(
    //         diaTargetFrames.first(),
    //         msReaderPointerAcc.ptr->getScanNumberVsScanTime(),
    //         PythiaParameterReader::genericPythiaParametersForTests(),
    //         diaTargetFrames.firstKey(),
    //         {10, 10},
    //         6,
    //         mzHashedVsCount,
    //         &msCalibratomatic,
    //         &turboXICMS1
    // );
    QCOMPARE(e, eNoError);

}

void CandidateScorertronTests::calculateScoresAndOtherStuffTooTest() {

    QSKIP("Fix tests for new API");

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
    QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    MsFrame msFrameMS1;
    TurboXIC turboXICMS1;

    e = buildInputData(
            &msReaderPointerAcc,
            &diaTargetFrames,
            &scanNumberVsScanPointsMS1,
            &ms1FramePtrs,
            &msFrameMS1,
            &turboXICMS1
    );
    QCOMPARE(e, eNoError);

    const QHash<MzHashed, int> mzHashedVsCount;
    MsCalibratomatic msCalibratomatic;

    const QString pythiaParamsFilePath
            = QDir(qApp->applicationDirPath()).filePath("test_params_wide_window.pythiaConfig");

    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
            pythiaParamsFilePath,
            &pythiaParameters
    );

    QMap<FrameIndex, ScanPoints> scanPoints;
    e = buildFrame(&scanPoints);
    QCOMPARE(e, eNoError);

    QMap<FrameIndex, ScanPoints*> scanPointsPtrs;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    for (auto it = scanPoints.begin(); it != scanPoints.end(); it++) {
        scanPointsPtrs.insert(it.key(), &it.value());
        scanNumberVsScanTime.insert(it.key(), it.key() / 100.0);
    }

    // CandidateScorertron candidateScorertron;
    // e = candidateScorertron.init(
    //         scanPointsPtrs,
    //         scanNumberVsScanTime,
    //         pythiaParameters,
    //         diaTargetFrames.firstKey(),
    //         {scanNumberVsScanTime.first(), scanNumberVsScanTime.last()},
    //         6,
    //         mzHashedVsCount,
    //         &msCalibratomatic,
    //         &turboXICMS1
    // );
    // QCOMPARE(e, eNoError);
    //
    // TargetDecoyCandidatePair targetDecoyCandidatePair = buildTargetDecoyCandidatePair();
    //
    // CandidateScores candidateScores;
    // e = candidateScorertron.calculateScores(
    //         targetDecoyCandidatePair.ms2IonsTarget(),
    //         &targetDecoyCandidatePair,
    //         &candidateScores
    //         );
    // QCOMPARE(e, eNoError);
    //
    // QCOMPARE(
    //         MathUtils::pRound(static_cast<double>(candidateScores.featuresArray[CandidateScores::Features::CosineSimSpectrumCubed]), 3),
    //         MathUtils::pRound(0.0767656341195, 3)
    //         );
    //
    // QHash<MzHashed , XICPoints> mzHashedVsXICPoints;
    // e = candidateScorertron.extractXICs(targetDecoyCandidatePair.ms2IonsTarget(), &mzHashedVsXICPoints);
    // QCOMPARE(mzHashedVsXICPoints.size(), 3);
    // QCOMPARE(mzHashedVsXICPoints.value(506335).size(), 6);
    // QCOMPARE(mzHashedVsXICPoints.value(407266).size(), 22);
    // QCOMPARE(mzHashedVsXICPoints.value(272125).size(), 31);

}


QTEST_MAIN(CandidateScorertronTests)

#include "CandidateScorertronTests.moc"
