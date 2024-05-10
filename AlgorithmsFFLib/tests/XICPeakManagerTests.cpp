//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>
#include <QtTest/QtTest>

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

        for (const MsScanInfo &msScanInfo : msScanInfos) {

            if (msScanInfo.msLevel < 2) {
                continue;
            }

            PythiaParameters pythiaParameters = PythiaParameterReader::genericPythiaParametersForTests();

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
    };

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

        XICPeakManager xicPeakManager;
        e = xicPeakManager.init(7, 2, 1, 0.5); ree;

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

        e = xicPeakManager.findPeaks(msFrame, mzValsToExtract, parallelInput.ppm); ree;

        if (parallelInput.cacheXICPeakManager) {
            const QString cacheFilePath = parallelInput.mzTargetKey + ".xicCache";
            e = xicPeakManager.cacheXICPeakManager(cacheFilePath); ree;
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
            = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML.prqFF");
            // = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML.prqFF");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath);
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(&diaTargetFrames);
    QCOMPARE(e, eNoError);


    const QString fragLibUri = "/home/anichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.fragLibFF";

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            fragLibUri,
            600,
            5000,
            &fragLibReaderRows
            );
    QCOMPARE(e, eNoError);

    fragLibReaderRows.resize(static_cast<int>(fragLibReaderRows.size() * 0.2));

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
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
        parallelInput.scanNumberVsScanTime = msReaderPointerAcc.ptr->getScanNumberVsScanTime();
        parallelInput.mzTargetKey = mzTargetKey;
        parallelInput.ppm = 20.0;
        parallelInput.topNFragmentIons = 6;
        parallelInput.cacheXICPeakManager = true;

        parallelInputs.push_back(parallelInput);
    }

    QFuture<Err> futures = QtConcurrent::mapped(parallelInputs, parallelLogic);
    futures.waitForFinished();

    for (const Err &res : futures) {
        e = res;
        QCOMPARE(e, eNoError);
    }

    QFuture<Err> futures2 = QtConcurrent::mapped(diaTargetFrames.keys(), readLogic);
    futures2.waitForFinished();
    for (const Err &res : futures2) {
        e = res;
        QCOMPARE(e, eNoError);
    }

    qDebug() << "Finsihed in" << et.elapsed() << "mSec";

}


QTEST_MAIN(XICPeakManagerTests)
#include "XICPeakManagerTests.moc"
