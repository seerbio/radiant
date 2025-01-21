#include <QtTest>
#include <QCoreApplication>

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "EigenUtils.h"
#include "MS2Ion.h"
#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
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

    // static TargetDecoyCandidatePair buildTargetDecoyCandidatePair() {
    //
    //     MS2Ion i1;
    //     i1.mz = 506.335;
    //     i1.intensity = 1.0;
    //     i1.charge = 1;
    //     i1.ionLabel = "y4";
    //
    //     MS2Ion i2;
    //     i2.mz = 407.266;
    //     i2.intensity = 0.488;
    //     i2.charge = 1;
    //     i2.ionLabel = "y3";
    //
    //     MS2Ion i3;
    //     i3.mz = 272.125;
    //     i3.intensity = 0.30;
    //     i3.charge = 1;
    //     i3.ionLabel = "b4";
    //
    //     QVector<MS2Ion> targetIons;
    //
    //     TargetDecoyCandidatePair targetDecoyCandidatePair(
    //             PeptideStringWithMods("DGVVLFK"),
    //             {i1, i2, i3},
    //             {i1, i2, i3},
    //             2,
    //             776.443,
    //             69.0608,
    //             3,
    //             0
    //             );
    //
    //     return targetDecoyCandidatePair;
    // }

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

            if (MathUtils::tZero(apexValue) || MathUtils::tZero(eVec.coeff(apexIndex))) {
                continue;
            }

            const float stopThreshold = apexValue * stopThresholdFraction;

            float rightStopVal = apexValue;
            int rightStopIndex = apexIndex;

            int rightCurrentIndex = apexIndex;
            while (rightCurrentIndex < eVec.size()) {

                const float currentValue = eVec.coeff(rightCurrentIndex);

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

    QSKIP("For testing");

    ERR_INIT

    const QString &intensityVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/intensity.csv");
    const QString &prodVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/prod.csv");
    const QString &limitsPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/limits.csv");

    QVector<float> prod;
    e = ObjectCSVWriters::readVectorFromFile(prodVecPath, &prod);
    QCOMPARE(e, eNoError);

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = simpleIntegrator(
        EigenUtils::convertQVectorToEigenVector(prod),
        0.65f,
        &peakIntegrationsVsIntensities
        );
    QCOMPARE(e, eNoError);

    QVector<QPair<int, int>> limits;
    for (const auto pr : peakIntegrationsVsIntensities) {
        limits.push_back(pr.first);
        qDebug() << pr << pr.first.first << pr.first.second;
    }

    e = ObjectCSVWriters::writeVectorToFile(limits, limitsPath);
    QCOMPARE(e, eNoError);

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
    //         MathUtils::pRound(static_cast<double>(candidateScores.featuresArray[Features::CosineSimSpectrumCubed]), 3),
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
