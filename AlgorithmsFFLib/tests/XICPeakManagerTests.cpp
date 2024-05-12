//
// Created by anichols on 11/07/2021.
//

#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>
#include <QtTest/QtTest>

#include "MsCalibratomatic.h"

class XICPeakManagerTests : public QObject
{
    Q_OBJECT

public:
    XICPeakManagerTests() = default;
    ~XICPeakManagerTests() override = default;

private Q_SLOTS:

    static void initTest();
    static void findPeaksTest();



private:


};


void XICPeakManagerTests::initTest() {

    ERR_INIT

    XICPeakManager xicPeakManager;

    e = xicPeakManager.init(3, 1, 1, 0.01);
    QCOMPARE(e, eNoError);

    e = xicPeakManager.init(2, 1, 1, 0.01);
    QCOMPARE(e, eError);

    e = xicPeakManager.init(3, 0, 1, 0.01);
    QCOMPARE(e, eError);

    e = xicPeakManager.init(3, 1, 0, 0.01);
    QCOMPARE(e, eError);

    e = xicPeakManager.init(3, 1, 1, 0.0);
    QCOMPARE(e, eError);

}

namespace {

    Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            const QVector<MsScanInfo> &msScanInfos,
            TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(msScanInfos); ree;
        e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;

        mzTargetKeyVsTargetDecoyCandidatePointers->clear();

        PythiaParameters pythiaParameters;
        e = PythiaParameterReader::buildPythiaParameters(
            "/home/anichols/Desktop/Data/ConfigFiles/test_params_decoys.pythiaConfig",
            &pythiaParameters
            );

        for (const MsScanInfo &msScanInfo : msScanInfos) {

            if (msScanInfo.msLevel < 2) {
                continue;
            }

            QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
            e = targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                    msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + pythiaParameters.precursorExtractionWindowThomsons),
                    msScanInfo.precursorTargetMz + (msScanInfo.isoWindowUpper + pythiaParameters.precursorExtractionWindowThomsons),
                    -1,
                    &targetDecoyPointers
                    ); ree;

            if (targetDecoyPointers.isEmpty()) {
                qWarning() << msScanInfo.targetKey() << "contains no target entries";
                continue;
            }

            mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyPointers);
            qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyPointers.size() << "targetDecoyPairs found";
        }

        ERR_RETURN
    }

    Err buildMzHashedVsCount(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
        int topNFragIons,
        QHash<MzHashed, int> *mzHashedVsCount
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyPointers); eee_absorb;

        for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {

            int counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsTarget()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;
            }

            counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsDecoy()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;
            }
        }

        ERR_RETURN
    }

    struct ParallelInput {
        QMap<ScanNumber, ScanPoints*> targetFrame;
        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairs;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        MzTargetKey mzTargetKey;
        float ppm = -1.0;
        int topNFragmentIons = -1;
        bool cacheXICPeakManager = false;
        float stopThresholdFraction = -1.0;
    };

    Err simpleIntegrator(
        const Eigen::VectorX<float> &vec,
        float stopThresholdFraction,
        QVector<QPair<PeakIntegrationIndexes, Intensity>> *peakIntegrationIndexesVsIntensity
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

                if (currentValue <= rightStopVal) {
                    rightStopVal = currentValue;
                    rightStopIndex = rightCurrentIndex;
                    rightCurrentIndex++;
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

                if (currentValue <= leftStopVal) {
                    leftStopVal = currentValue;
                    leftStopIndex = leftCurrentIndex;
                    leftCurrentIndex--;
                    continue;
                }

                break;
            }

            peakIntegrationIndexesVsIntensity->push_back({
                 {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size() - 1))},
                 apexValue}
            );

            for (int i = leftStopIndex; i <= rightStopIndex; i++) {
                eVec.coeffRef(i) = 0.0;
            }
        }

        std::sort(
            peakIntegrationIndexesVsIntensity->rbegin(),
            peakIntegrationIndexesVsIntensity->rend(),
            [](const QPair<PeakIntegrationIndexes, Intensity> &l, const QPair<PeakIntegrationIndexes, Intensity> &r) {
                return l.second < r.second;
            }
            );

        ERR_RETURN
    }

    Err buildMatriciesAndIntegrationVector(
        const QVector<MS2Ion> &ms2Ions,
        const Eigen::VectorX<float> &kernel,
        const Eigen::VectorX<float> &kernelIntegration,
        XICPeakManager *xicPeakManager,
        Eigen::MatrixX<float> *matOutput,
        Eigen::VectorX<float> *integrationVec
        ) {

        ERR_INIT

        FrameIndex frameIndexMax = -1;

        QVector<XICPoints> xicPointses;
        xicPointses.reserve(ms2Ions.size());
        for (const MS2Ion &ms2Ion : ms2Ions) {
            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

            if (xicPoints.empty()) {
                xicPointses.push_back(XICPoints());
                continue;
            }

            xicPointses.push_back(xicPoints);

            FrameIndex frameIndexMaxMs2Ion = std::max_element(
                xicPoints.begin(),
                xicPoints.end(),
                [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
                )->scanNumber;

            frameIndexMax = std::max(frameIndexMax, frameIndexMaxMs2Ion);
        }

        Eigen::MatrixX<float> mat(frameIndexMax + 1, xicPointses.size());
        mat.setZero();
        for (int col = 0; col < xicPointses.size(); col++) {
            const XICPoints &xicPoints = xicPointses.at(col);
            for (const XICPoint &xicPoint : xicPoints) {
                mat.coeffRef(xicPoint.scanNumber, col) = xicPoint.intensity;
            }
        }
        mat = EigenKernelUtils::applyKernelToEachColumnInMatrix(mat, kernel);

        Eigen::MatrixX<float> matCount = mat;
        matCount = (matCount.array() > 0.1).select(1.0, matCount);
        EigenUtils::thresholdMatrix(0.0f, &matCount);

        Eigen::VectorX<float> integrationVecLocal = matCount.rowwise().sum();
        constexpr float minPeakCount = 3.9;
        EigenUtils::thresholdVector(minPeakCount, &integrationVecLocal);

        EigenUtils::thresholdMatrix(0.0f, &mat);
        *matOutput = mat;
        *integrationVec = EigenKernelUtils::convolveVectorWithKernel(integrationVecLocal, kernelIntegration);

        ERR_RETURN
    }

    namespace {

        QVector<int> findStartApexes(
            const Eigen::MatrixX<float> &matBlockApexes,
            int apexIndex
            ) {

            QVector<int> apexesStartsToUse;
            apexesStartsToUse.reserve(static_cast<int>(matBlockApexes.cols()));
            for (int col = 0; col < matBlockApexes.cols(); col++) {

                const Eigen::VectorX<float> &apexColumn = matBlockApexes.col(col);

                if (apexColumn.coeff(apexIndex) > 0) {
                    apexesStartsToUse.push_back(apexIndex);
                    continue;
                }

                for (int rowFromCenter = 1; rowFromCenter < apexColumn.size(); rowFromCenter++) {

                    const int rowLeftIndex = apexIndex - rowFromCenter;
                    const int rowRightIndex = apexIndex + rowFromCenter;

                    int rowLeftIndexValue = -1;
                    int rowRightIndexValue = -1;

                    if (rowLeftIndex >= 0) {
                        rowLeftIndexValue = static_cast<int>(apexColumn.coeff(rowLeftIndex));
                    }

                    if (rowRightIndex < apexColumn.size()) {
                        rowRightIndexValue = static_cast<int>(apexColumn.coeff(rowRightIndex));
                    }

                    if (rowLeftIndexValue > 0 && rowRightIndexValue > 0) {
                        const int higherIndex = rowLeftIndexValue >= rowRightIndexValue ? rowLeftIndex : rowRightIndex;
                        apexesStartsToUse.push_back(higherIndex);
                        break;
                    }
                    else if (rowLeftIndexValue > 0) {
                        apexesStartsToUse.push_back(rowLeftIndex);
                        break;
                    }
                    else if (rowRightIndexValue > 0) {
                        apexesStartsToUse.push_back(rowRightIndex);
                        break;
                    }

                    if (rowLeftIndex <= 0 && rowRightIndex >= apexColumn.size() - 1) {
                        apexesStartsToUse.push_back(-1);
                        break;
                    }
                }
            }

            return apexesStartsToUse;
        }

        QPair<int, int> simpleIntegratorTrimmer(
            const Eigen::VectorX<float> &vec,
            int apexRowIndex,
            float stopThresholdFraction
            ) {

            const float nearZero = 0.001;
            const float apexIndexValue = vec.coeff(apexRowIndex);
            const float stopThresholdValue = apexIndexValue * stopThresholdFraction;

            int rightStopIndex = apexRowIndex;
            int rightCurrentIndex = apexRowIndex;
            while (rightCurrentIndex < vec.size()) {

                const float currentValue = vec.coeff(rightCurrentIndex);
                if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                    rightStopIndex = rightCurrentIndex;
                    break;
                }

                if (currentValue > nearZero) {
                    rightStopIndex = rightCurrentIndex;
                    rightCurrentIndex++;
                    continue;
                }

                break;
            }

            int leftStopIndex = apexRowIndex;
            int leftCurrentIndex = apexRowIndex;
            while (leftCurrentIndex >= 0) {

                const float currentValue = vec(leftCurrentIndex);
                if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                    leftStopIndex = leftCurrentIndex;
                    break;
                }

                if (currentValue > nearZero) {
                    leftStopIndex = leftCurrentIndex;
                    leftCurrentIndex--;
                    continue;
                }

                break;
            }

            return {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size()) - 1)};
        }

    }//namespace
    Err parallelLogic (const ParallelInput &parallelInput) {

        ERR_INIT

            QElapsedTimer et;
            et.start();

            QHash<MzHashed, int> mzHashedVsCount;
            e = buildMzHashedVsCount(
                parallelInput.targetDecoyCandidatePairs,
                parallelInput.topNFragmentIons,
                &mzHashedVsCount
                ); ree;


            MsFrame msFrame;
            e = msFrame.init(parallelInput.targetFrame, parallelInput.scanNumberVsScanTime); ree;

            const QList<int> &mzHashedVsCountKeys = mzHashedVsCount.keys();
            QVector<float> mzValsToExtract;
            std::transform(
                mzHashedVsCountKeys.begin(),
                mzHashedVsCountKeys.end(),
                std::back_inserter(mzValsToExtract),
                [](int mzHashed){return MathUtils::unHashDecimal<float>(mzHashed, S_GLOBAL_SETTINGS.HASHING_PRECISION);}
                );

            XICPeakManager xicPeakManager;
            e = xicPeakManager.init(msFrame, mzValsToExtract, parallelInput.ppm); ree;

            if (parallelInput.cacheXICPeakManager) {
                const QString cacheFilePath = parallelInput.mzTargetKey + ".xicCache";
                e = xicPeakManager.cacheXICPeakManager(cacheFilePath); ree;
            }

            Eigen::MatrixX<float> kernel;
            e = EigenKernelUtils::buildSavitzkyGolayKernel(5, 2, 0, 1, &kernel); ree;
            Eigen::VectorX<float> kernelVec(kernel);

            Eigen::MatrixX<float> kernelIntegration;
            e = EigenKernelUtils::buildSavitzkyGolayKernel(7, 2, 0, 1, &kernelIntegration); ree;
            Eigen::VectorX<float> kernelIntegrationVec(kernelIntegration);

            for (TargetDecoyCandidatePair *tdcp : parallelInput.targetDecoyCandidatePairs) {

                QVector<MS2Ion> ms2IonsTarget = tdcp->ms2IonsTarget();
                ms2IonsTarget.resize(std::min(parallelInput.topNFragmentIons, ms2IonsTarget.size()));

                Eigen::VectorX<float> integrationVec;
                Eigen::MatrixX<float> mat;
                e = buildMatriciesAndIntegrationVector(
                    ms2IonsTarget,
                    kernelVec,
                    kernelIntegration,
                    &xicPeakManager,
                    &mat,
                    &integrationVec
                    ); ree;

                QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationIndexesVsIntensity;
                e = simpleIntegrator(
                    integrationVec,
                    parallelInput.stopThresholdFraction,
                    &peakIntegrationIndexesVsIntensity
                    ); ree;

                peakIntegrationIndexesVsIntensity.resize(std::min(10, peakIntegrationIndexesVsIntensity.size()));
                for (const QPair<PeakIntegrationIndexes, Intensity> &pr : peakIntegrationIndexesVsIntensity) {

                    Eigen::MatrixX<float> matBlock = mat.block(
                        pr.first.first,
                        0,
                        pr.first.second - pr.first.first + 1,
                        mat.cols()
                        ).eval();

                    if (pr.first.second == pr.first.first) {
                        continue;
                    }

                    QVector<QVector<int>> apexIndexesByColumn;
                    apexIndexesByColumn.reserve(matBlock.cols());
                    for (int col = 0; col < matBlock.cols(); col++) {
                        const Eigen::VectorX<float> &column = matBlock.col(col);
                        apexIndexesByColumn.push_back(EigenUtils::apexesIndexesOnly(column));
                    }

                    Eigen::MatrixX<float> matBlockApexes(matBlock.rows(), matBlock.cols());
                    matBlockApexes.setZero();
                    for (int col = 0; col < apexIndexesByColumn.size(); col++) {
                        const QVector<int> &columnApexes = apexIndexesByColumn.at(col);
                        for (int row : columnApexes) {
                            matBlockApexes.coeffRef(row, col) = matBlock.coeff(row, col);
                        }
                    }

                    const Eigen::VectorX<float> integrationVecSegment = integrationVec.segment(pr.first.first, pr.first.second - pr.first.first + 1);

                    const QPair<int, float> apexIndex = EigenUtils::returnTopIndexAndValue(integrationVecSegment);

                    // if (pr.first.first > 0) {
                    //     continue;
                    // }
                    //
                    // // std::cout << integrationVec << std::endl;
                    // std::cout << peakIntegrationIndexesVsIntensity.size() << std::endl;
                    // std::cout << apexIndex.first << " " << apexIndex.second << std::endl;
                    // std::cout << Eigen::RowVectorX<float>(integrationVecSegment) << std::endl;


                    const QVector<int> apexStarts = findStartApexes(matBlockApexes, apexIndex.first);
                    e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockApexes.cols()));
                    if (e != eNoError) {
                        std::cout << apexIndex.first << " " << apexIndex.second << " " << pr.first.first << " " << pr.first.second << std::endl;;
                        std::cout << Eigen::RowVectorX<float>(integrationVecSegment) << std::endl;
                        std::cout << "-------" << std::endl;
                        std::cout << matBlockApexes << std::endl;
                        std::cout << "***********" << std::endl;
                    }
                    Eigen::MatrixX<float> matBlockTrimmedColumns(matBlock.rows(), matBlock.cols());
                    matBlockTrimmedColumns.setZero();
                    for (int col = 0; col < matBlock.cols(); col++) {

                        const Eigen::VectorX<float> &colVec = matBlock.col(col);
                        const QPair<int, int> peakLimits = simpleIntegratorTrimmer(
                            colVec,
                            apexStarts.at(col),
                            parallelInput.stopThresholdFraction
                            );

                        for(int row = peakLimits.first; row <= peakLimits.second; row++) {
                            matBlockTrimmedColumns.coeffRef(row, col) = colVec.coeff(row);
                        }
                    }

                    for (int a : apexStarts) {
                        std::cout << a << std::endl;
                    }
                    std::cout << apexIndex.first << " " << apexIndex.second << std::endl;;
                    std::cout << Eigen::RowVectorX<float>(integrationVecSegment) << std::endl;
                    std::cout << "-------" << std::endl;
                    std::cout << matBlock << std::endl;
                    std::cout << "-------" << std::endl;
                    std::cout << matBlockApexes << std::endl;
                    std::cout << "-------" << std::endl;
                    std::cout << matBlockTrimmedColumns << std::endl;
                    std::cout << "***********" << std::endl;

                }

            }

            qDebug() << parallelInput.mzTargetKey << et.elapsed();

        ERR_RETURN
    }

    Err readLogic(const MzTargetKey &mzTargetKey) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        const QString cacheFilePath = mzTargetKey + ".xicCache";
        XICPeakManager xicPeakManager;
        e = xicPeakManager.loadXICPeakManagerCache(cacheFilePath);

        qDebug() << et.elapsed() << mzTargetKey;

        ERR_RETURN
    }

}//namespace
void XICPeakManagerTests::findPeaksTest() {

    ERR_INIT

    const QString msDataFilePath
            // = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML.prqFF");
            // = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML.prqFF");
            = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML.prqFF");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(&diaTargetFrames);
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, ScanTime> scanNumberVsScanTime = msReaderPointerAcc.ptr->getScanNumberVsScanTime();

    const QString fragLibUri
            = "/home/anichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF";
            // = "/home/anichols/Desktop/Data/Libraries/predicted_lib.speclib.fragLibFF";

    PythiaParameters pythiaParameters;
    e = PythiaParameterReader::buildPythiaParameters(
        "/home/anichols/Desktop/Data/ConfigFiles/test_params_decoys.pythiaConfig",
        &pythiaParameters
        );
    QCOMPARE(e, eNoError);

    const double massMin
        = pythiaParameters.peptideLengthMin * Molecule(MolecularFormulas::alanineFormula).monoisotopicMass();

    const double massMax
            = pythiaParameters.peptideLengthMax * Molecule(MolecularFormulas::tryptophanFormula).monoisotopicMass();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            fragLibUri,
            massMin,
            massMax,
            &fragLibReaderRows
            );
    QCOMPARE(e, eNoError);

    // fragLibReaderRows.resize(static_cast<int>(fragLibReaderRows.size() * 0.2));
    fragLibReaderRows.resize(static_cast<int>(1e4));
    // fragLibReaderRows.resize(pythiaParameters.trancheSizeMax);

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            pythiaParameters,
            &fragLibReaderRows
            );
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        msReaderPointerAcc.ptr->getUniqueTandemMsScanInfos(),
        &targetDecoyCandidatePairManager,
        &mzTargetKeyVsTargetDecoyCandidatePointers
        );
    QCOMPARE(e, eNoError);

    QElapsedTimer et;
    et.start();

    QVector<ParallelInput> parallelInputs;
    for (auto it = mzTargetKeyVsTargetDecoyCandidatePointers.begin(); it != mzTargetKeyVsTargetDecoyCandidatePointers.end(); ++it) {

        const MzTargetKey &mzTargetKey = it.key();
        ParallelInput parallelInput;
        parallelInput.targetFrame = diaTargetFrames.value(mzTargetKey);
        parallelInput.targetDecoyCandidatePairs = it.value();
        parallelInput.scanNumberVsScanTime = scanNumberVsScanTime;
        parallelInput.mzTargetKey = mzTargetKey;
        parallelInput.ppm = 20.0;
        parallelInput.topNFragmentIons = 6;
        parallelInput.cacheXICPeakManager = false;
        parallelInput.stopThresholdFraction = 0.5;

        parallelInputs.push_back(parallelInput);
    }

// #define RUNP
#ifdef RUNP
    et.restart();
    QFuture<Err> futures = QtConcurrent::mapped(parallelInputs, parallelLogic);
    futures.waitForFinished();

    for (const Err &res : futures) {
        e = res;
        QCOMPARE(e, eNoError);
    }
#else
    for (const ParallelInput &pi : parallelInputs) {
        parallelLogic(pi);
    }
#endif

    // qDebug() << et.restart() << "write";

    // et.restart();
    // QFuture<Err> futures2 = QtConcurrent::mapped(diaTargetFrames.keys(), readLogic);
    // futures2.waitForFinished();
    // for (const Err &res : futures2) {
    //     e = res;
    //     QCOMPARE(e, eNoError);
    // }
    //
    // qDebug() << et.restart() << "read";

    qDebug() << "Finsihed in" << et.elapsed() << "mSec";
}


QTEST_MAIN(XICPeakManagerTests)
#include "XICPeakManagerTests.moc"
