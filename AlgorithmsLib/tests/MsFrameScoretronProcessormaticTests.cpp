#include "MsFrameScoretronProcessormatic.h"

#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsFrameScoretron.h"
#include "FrameExtractsReader.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"

#include <QElapsedTimer>
#include <QVector>
#include <QtTest>


class MsFrameScoretronProcessormaticTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronProcessormaticTests() = default;
    ~MsFrameScoretronProcessormaticTests() override = default;


private slots:

    void rescoreMsFrameTest();


};


void MsFrameScoretronProcessormaticTests::rescoreMsFrameTest() {

    ERR_INIT

    const QString frameExtractionFilepath
        = "/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.444952.frameExtractions";

    QVector<FrameExtractsReaderRow> frameExtracts;
    e = ParquetReader::read(
            frameExtractionFilepath,
            &frameExtracts
            );
    QCOMPARE(e, eNoError);

    MsFrameScoretronProcessormatic msFrameScoretronProcessormatic;
    e = msFrameScoretronProcessormatic.init(PythiaParameterReader::genericPythiaParametersForTests());
    QCOMPARE(e, eNoError);


    std::sort(
            frameExtracts.begin(),
            frameExtracts.end(),
            [](const FrameExtractsReaderRow &l, const FrameExtractsReaderRow &r){return l.frameIndexApex < r.frameIndexApex;}
    );

    QMap<int, bool> debugged;
    QVector<FrameExtractsReaderRow> frameIndexFrameExtractsReaderRows;
    for (const FrameExtractsReaderRow &fer : frameExtracts) {

        if (!debugged.value(fer.scanNumberApex)) {
            qDebug() << fer.scanNumberApex << fer.frameIndexApex;
        }
        debugged.insert(fer.scanNumberApex, true);

        if (fer.scanNumberApex != 8261) {
            continue;
        }

        qDebug() << fer.peptideStringWithMods;
        frameIndexFrameExtractsReaderRows.push_back(fer);
    }


    qDebug() << "candidate size" << frameIndexFrameExtractsReaderRows.size();

    e = msFrameScoretronProcessormatic.deconvolveScanCandidate(frameIndexFrameExtractsReaderRows);
    QCOMPARE(e, eNoError);

//        ScanPoints scanPoints;
//        e = ParallelUtils::zip(
//                fer.yIonMzValsActual,
//                fer.yIonIntesitiesActual,
//                &scanPoints
//                );
//        QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsFrameScoretronProcessormaticTests)

#include "MsFrameScoretronProcessormaticTests.moc"
